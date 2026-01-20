module S5Grid

using Base, Dates, Printf
using NCDatasets
using Statistics
using Glob
using JSON
using ProgressMeter

export run_gridding, run_gridding_km

### some helper functions...

# Calculate number of bins for a specific latitude row to maintain km resolution
function get_n_lon(lat_val, target_km)
    # Earth circumference at equator approx 40075 km
    # At this latitude:
    circ_at_lat = 40075.0 * cos(deg2rad(lat_val))
    
    # How many bins of target_km fit?
    n_bins = ceil(Int, circ_at_lat / target_km)
    return max(1, n_bins) # Ensure at least 1 bin
end

function getFilter(name, jsonDict)
    ff = try
        jsonDict[name]
    catch
        Dict()
    end
    return ff
end

# This splits up the entire region into one grid set
function divLine!(lat1,lon1,lat2,lon2,n, points,j )
    dLat     = (lat2-lat1)/(2n)
    dLon     = (lon2-lon1)/(2n)
    startLat = lat1+dLat
    startLon = lon1+dLon
    @inbounds for i in 1:n
        #println(startLat+2*(i-1)*dLat, " ", startLon+2*(i-1)*dLon)
        #weights[(iLat-minLat+1), (iLon-minLon+1)]+=1
        points[j,i,1] = startLat+2*(i-1)*dLat
        points[j,i,2] = startLon+2*(i-1)*dLon
    end
end

# For first run (2 baselines)
function divLine2!(lat1,lon1,lat2,lon2,n, lats, lons)
    dLat     = (lat2-lat1)/(2n)
    dLon     = (lon2-lon1)/(2n)
    startLat = lat1+dLat
    startLon = lon1+dLon
    @inbounds for i in 1:n
        lats[i] = startLat+2*(i-1)*dLat
        lons[i] = startLon+2*(i-1)*dLon
    end
end

# Divide each polygon into multiple points
function getPoints!(points, vert_lat, vert_lon, n,lats_0, lons_0,lats_1, lons_1 )
    # Get reference points for two lines at the extremes:
    # println(vert_lat)
    divLine2!(vert_lat[1],vert_lon[1],vert_lat[2],vert_lon[2],n,lats_0, lons_0)
    divLine2!(vert_lat[4],vert_lon[4],vert_lat[3],vert_lon[3],n,lats_1, lons_1)
    for i in 1:n
         divLine!(lats_0[i], lons_0[i] ,lats_1[i], lons_1[i],n, points,i )
    end
end

## this fixes a bug when boundaries coordinates are at root level:
function getNC_var(fin, path, DD::Bool)
    loc = split(path, r"/")
    # 1. Get the variable object regardless of depth
    var_obj = nothing
    if length(loc) == 1
        var_obj = fin[path]
    else
        gr = fin
        for i in 1:length(loc)-1
            gr = gr.group[loc[i]]
        end
        var_obj = gr[loc[end]]
    end

    # 2. Unified Reshape Logic
    si = size(var_obj)
    raw_data = var_obj.var[:]

    if DD
        # Handle 4-corner boundary pixels
        if si[1] == 4
            return reshape(raw_data, 4, prod(si[2:end]))'
        elseif si[end] == 4
            return reshape(raw_data, prod(si[1:end-1]), 4)
        else
            return reshape(raw_data, prod(si[1:end-1]), si[end]) # Fallback
        end
    else
        # Handle standard 2D/1D data
        return reshape(raw_data, prod(si))
    end
end


function getNC_attrib(fin, path, attri)
    loc = split(path ,r"/")
    #println(loc)
    if length(loc)==1
        return fin[path].attrib[attri]
    elseif length(loc)>1
        gr = []
        for i in 1:length(loc)-1
            if i==1
                gr = fin.group[loc[i]]
            else
                gr = gr.group[loc[i]]
            end
        end
        return gr[loc[end]].attrib[attri]
    end
end



function floorPoints!(points,ix,iy,n)
    pointsX = @view points[:,:,1]
    pointsY = @view points[:,:,2]
    for iC in LinearIndices(pointsX)
        ix[iC]=floor(pointsX[iC])
        iy[iC]=floor(pointsY[iC])
    end
end


function favg_all!(arr,std_arr, weight_arr, compSTD, lat,lon,inp,s,s2,n, latMin, latMax, lonMin, lonMax, nLat, nLon, points)

    dim = size(arr)
    #println(dim)
    # Predefine some arrays to reduce allocations
    ix      = zeros(Int32,n^2)
    iy      = zeros(Int32,n^2)
    lats_0  = zeros(n)
    lons_0  = zeros(n)
    lats_1  = zeros(n)
    lons_1  = zeros(n)
    iLon    = floor.(Int32,lon)
    iLat    = floor.(Int32,lat)
    minLat  = minimum(Int32,floor.(lat), dims=2)
    maxLat  = maximum(Int32,floor.(lat), dims=2)
    minLon  = minimum(Int32,floor.(lon), dims=2)
    maxLon  = maximum(Int32,floor.(lon), dims=2)
    distLon = maxLon-minLon
    # How many individual grid cells might actually be there:
    dimLat = maxLat-minLat
    dimLon = maxLon-minLon
    fac = Float32(1/n^2)
    #@show s
    
    
    @inbounds for i in 1:s
        #println(i, " ", dimLat[i], " ", dimLon[i])
        # Take it easy if all corners already fall into one grid box:
        #@show distLon[i]
        if (dimLat[i]==0) & (dimLon[i]==0)
            
            weight_arr[iLon[i,1],iLat[i,1]] += 1
            for z in 1:s2
                mean_old = arr[iLon[i,1],iLat[i,1],z];
                arr[iLon[i,1],iLat[i,1],z] = mean_old .+ 1/weight_arr[iLon[i,1],iLat[i,1]] .* (inp[i,z]-mean_old);
                if compSTD
                    std_arr[iLon[i,1],iLat[i,1],z] += (inp[i,z]-mean_old) .* (inp[i,z]-arr[iLon[i,1],iLat[i,1],z])
                end
            end
        # if not, compute appropriate weights
        elseif (distLon[i])<n
            getPoints!(points,lat[i,:],lon[i,:],n,lats_0, lons_0,lats_1, lons_1 )

            floorPoints!(points,ix,iy,n)
            #ix .= floor.(Int32,points[:,:,1][:])
            #iy .= floor.(Int32,points[:,:,2][:])
            
            @inbounds for j in eachindex(ix)
                weight_arr[iy[j],ix[j]] += fac;
                for z in 1:s2
                    mean_old = arr[iy[j],ix[j],z]
                    arr[iy[j],ix[j],z] = mean_old + fac/weight_arr[iy[j],ix[j]] * (inp[i,z]-mean_old);
                    if compSTD
                        std_arr[iy[j],ix[j],z] +=  fac * (inp[i,z]-mean_old) * (inp[i,z]-arr[iy[j],ix[j],z])
                    end
                end
            end
        end
    end
end


# Helper to handle the coordinate transformation for a single point
# Returns (0,0) if out of bounds
@inline function get_km_indices(lat_val, lon_val, latMin, dLat, nLat, n_lons_per_row, lonMin, lonMax)
    # 1. Calculate Row (Latitude)
    # Note: 1-based indexing for Julia
    row_idx = floor(Int, (lat_val - latMin) / dLat) + 1
    
    if row_idx < 1 || row_idx > nLat
        return (0, 0)
    end

    # 2. Calculate Column (Longitude) based on THIS row's resolution
    n_bins = n_lons_per_row[row_idx]
    
    # Normalize longitude to 0..1 range
    # Ensure lon_val is within [lonMin, lonMax] or handle wrapping
    lon_span = lonMax - lonMin
    norm_lon = (lon_val - lonMin) / lon_span

    # Map to bin index
    col_idx = floor(Int, norm_lon * n_bins) + 1
    
    # Handle wrap-around or edge precision issues
    if col_idx > n_bins 
        col_idx = 1 # Wrap 180 -> -180 if needed, or clamp
    elseif col_idx < 1
        col_idx = n_bins 
    end

    return (row_idx, col_idx)
end

function favg_all_km!(arr, std_arr, weight_arr, compSTD, 
                      lat, lon, inp, s, s2, nGrid, 
                      latMin, latMax, lonMin, lonMax, 
                      nLat, n_lons_per_row, points)

    # Calculate dLat once (assumed constant in degrees)
    dLat = (latMax - latMin) / nLat
    
    # Sub-pixel weighting factor
    fac = Float32(1.0 / nGrid^2)
    
    # Pre-allocate scratch arrays for the polygon corners (from original code)
    lats_0 = zeros(nGrid)
    lons_0 = zeros(nGrid)
    lats_1 = zeros(nGrid)
    lons_1 = zeros(nGrid)

    @inbounds for i in 1:s
        # Get bounding box of the satellite pixel
        # Note: Original code logic for min/max/dim is simplified here for clarity
        # assuming 'lat' and 'lon' inputs to this function are the CENTER points 
        # and we use getPoints! to find the spread.
        
        # If you have bounding box info passed in, use it to decide if splitting is needed.
        # always split if nGrid > 1
        # --- CASE A: Simple Center Point (Fast Mode or Pixel is tiny) ---
        # NOTE: nGrid is currently hardcoded 
        if nGrid == 1
            r, c = get_km_indices(lat[i], lon[i], latMin, dLat, nLat, n_lons_per_row, lonMin, lonMax)
            
            if r > 0 # valid index
                weight_arr[r][c] += 1
                for z in 1:s2
                    mean_old = arr[r][c, z]
                    diff = inp[i, z] - mean_old
                    arr[r][c, z] = mean_old + (diff / weight_arr[r][c])
                    
                    if compSTD
                        std_arr[r][c, z] += diff * (inp[i, z] - arr[r][c, z])
                    end
                end
            end

        # --- CASE B: Sub-pixel Splitting (High Accuracy) ---
        else
            # Assuming the calling function passed the corner coordinates 
            # e.g., lat is actually a matrix [pixels, 4] of corners
            getPoints!(points, lat[i,:], lon[i,:], nGrid, lats_0, lons_0, lats_1, lons_1)
            
            # Now iterate over the generated sub-points
            # Note: cannot use floorPoints! anymore because it assumes rectangular grid.
            # Instead calculate indices for every sub-point.
            pointsX = @view points[:,:,1] # Lats
            pointsY = @view points[:,:,2] # Lons
            
            for sub_p in eachindex(pointsX)
                # Map this specific sub-point to the Ragged Grid
                pLat = pointsX[sub_p]
                pLon = pointsY[sub_p]
                
                r, c = get_km_indices(pLat, pLon, latMin, dLat, nLat, n_lons_per_row, lonMin, lonMax)
                
                if r > 0 # If valid
                    # Add fractional weight
                    weight_arr[r][c] += fac
                    
                    # Welford/Mean Update
                    # Note: Using 'fac' as weight in the update is slightly different 
                    # than standard Welford but standard for "splatting" grids.
                    for z in 1:s2
                        mean_old = arr[r][c, z]
                        # weight the update by the accumulation so far
                        curr_w = weight_arr[r][c]
                        
                        # Update Mean
                        arr[r][c, z] = mean_old + (fac / curr_w) * (inp[i, z] - mean_old)
                        
                        if compSTD
                            # Variance update is tricky with fractional weights, 
                            # simpler approximation:
                            std_arr[r][c, z] += fac * (inp[i, z] - mean_old) * (inp[i, z] - arr[r][c, z])
                        end
                    end
                end
            end
        end
    end
end


# Convert 0–360° to -180–180° if needed
fixlon(x) = ifelse.(x .> 180, x .- 360, x)


function run_gridding(files_any::Any;
        json_dict::String = "jsonFiles/S5_O3_dict.json",
        monthly::Bool=false, 
        compSTD::Bool=false, 
        latMin::Float64=-90.0, 
        latMax::Float64=90.0, 
        lonMin::Float64=-180.0, 
        lonMax::Float64=180.0,
        dLat::Float64=1.0, 
        dLon::Float64=1.0,
        startDate_str::String="2025-11-25", 
        stopDate_str::String="2015-11-26",
        dDays::Int64=8, 
        oversample_temporal::Float64=1.0,
        FillValue=-999
    )

    println("Running S5 gridding (Python mode, in-memory, rectangular grid)")
    
    # Convert PyList{Any} or Any other list to a Julia Vector{String}
    # at the same time convert list into a Set for fast lookup 
    files_set = Set(convert(Vector{String}, collect(files_any)))

    startDate = DateTime(startDate_str)
    stopDate = DateTime(stopDate_str)
    dDay  = monthly ? Dates.Month(dDays) : Dates.Day(dDays)
    println(startDate, " ", stopDate)
    cT = length(startDate:dDay:stopDate)
    
    eps = dLat/100

    lat = collect(latMin+dLat/2.:dLat:latMax-dLat/2.0+eps)
    lon = collect(lonMin+dLon/2.:dLon:lonMax-dLon/2.0+eps)
    println("Output dimension (time/lon/lat):")
    println(cT, "/", length(lon),"/", length(lat))

    # Define gridded variables:
    n=zeros(Float32,(length(lat),length(lon)))
    
    # Parse JSON files as dictionary
    jsonDict = JSON.parsefile(json_dict)
    d2       = jsonDict["basic"]
    dGrid    = jsonDict["grid"]

    # Read all filters:
    f_eq = getFilter("filter_eq",jsonDict)
    f_gt = getFilter("filter_gt",jsonDict)
    f_lt = getFilter("filter_lt",jsonDict)

    # Get file naming pattern (needs YYYY MM and DD in there)
    fPattern = jsonDict["filePattern"]
    # Get main folder for files:
    # folder   = jsonDict["folder"] # not needed as files are given via python... files already exist
    
    # hard-coded value to divide each polygon, can be changed:
    nGrid = 10;
    #println("nGrid is set to: ", nGrid)
    #####

    points = zeros(Float32,(nGrid,nGrid,2))

    # Loop through time:
    # Time counter
    cT = 1
    p1 = Progress(cT)
    mat_data          = zeros(Float32,(length(lon),length(lat),length(dGrid)))
    mat_data_variance = zeros(Float32,(length(lon),length(lat),length(dGrid)))
    mat_data_weights  = zeros(Float32,(length(lon),length(lat)))
    nTime = length(startDate:dDay:stopDate)

    dates = collect(startDate:dDay:stopDate)

    println("Creating output arrays...")
    out_mat_av = zeros(Float32,(length(dates),length(lon),length(lat),length(dGrid)))
    out_mat_n  = zeros(Float32,(length(dates),length(lon),length(lat)))
    out_mat_std = compSTD ? zeros(Float32,(length(dates),length(lon),length(lat),length(dGrid))) : nothing # only needed if compSTD is true

    for cT in eachindex(dates)
        d = dates[cT]
    
        println("Gridding time slice ", d, " (",cT,"/", nTime,")")
        ProgressMeter.next!(p1; showvalues = [(:Time, d)])
    
        # 1. Define the search window
        # This matches the range from orgiginal code: [d, d + dDay*oversample_temporal - 1 Day]
        date_range = d:Dates.Day(1):d + dDay * oversample_temporal - Dates.Day(1)
        
        ### find files within this time step from provided list
        # Create specific patterns for each date
        target_patterns = [replace(fPattern, 
                            "YYYY" => Dates.format(di, "yyyy"), 
                            "MM"   => Dates.format(di, "mm"), 
                            "DD"   => Dates.format(di, "dd")) 
                           for di in date_range]
        
        # Create a single Regex that combines all dates: e.g., (pattern1|pattern2|pattern3)
        # replace '*' with '.*' to make it a valid regex wildcard
        combined_regex = Regex(join(replace.(target_patterns, "*" => ".*"), "|"))
        
        # Filter once
        files_tmp = filter(f -> occursin(combined_regex, f), files_set)
        # Convert the Set to a Vector (Array) so it supports positional and boolean indexing.
        files_vector = collect(files_tmp)        
        
        println("  -> Files found (", length(files_vector), "): ")#, files)
        #######

        fileSize = Int[];
        for f in files_vector # Iterate over the Vector
            fileSize = [fileSize;stat(f).size]
        end

        println("File size: ", fileSize) 
        
        # Loop through all files within this time step
        p = Progress(length(files))

        for a in files_vector[fileSize.>0]
            try
                fin = Dataset(a)
                println("Read, ", a)
                
                lat_center = getNC_var(fin, d2["lat"],false)
                lon_center = getNC_var(fin, d2["lon"],false)

                # --- convert 0–360 to –180–180 (center coords) if needed ---
                lon_center = fixlon(lon_center)
    
                bool_found = (lat_center.>latMin) .+ (lat_center.<latMax) .+ (lon_center.>lonMin) .+ (lon_center.<lonMax)
                if any(bool_found.==4)
    
                    # Read lat/lon bounds (required, maybe can change this to simple gridding in the future with just center):
                    lat_in_ = getNC_var(fin, d2["lat_bnd"],true)
                    lon_in_ = getNC_var(fin, d2["lon_bnd"],true)
                    
                    #println("Read")
                    dim = size(lat_in_)
    
                    # Transpose if orders are swapped
                    if dim[1]==4
                        lat_in_ = lat_in_'
                        lon_in_ = lon_in_'
                    end
    
                    # --- convert 0–360 to –180–180 (bounds) if needed ---
                    lon_in_ = fixlon(lon_in_)

                    # Find all indices within lat/lon bounds:
                    minLat = minimum(lat_in_, dims=2)
                    maxLat = maximum(lat_in_, dims=2)
                    minLon = minimum(lon_in_, dims=2)
                    maxLon = maximum(lon_in_, dims=2)
    
                    # Get indices within the lat/lon bounding box and check filter criteria (the last one filters out data crossing the date boundary):
                    bool_add = (minLat[:,1].>latMin) .+ (maxLat[:,1].<latMax) .+ (minLon[:,1].>lonMin) .+ (maxLon[:,1].<lonMax) .+ ((maxLon[:,1].-minLon[:,1]).<50)
                    
                    bCounter = 5
                    # Look for equalities
                    for (key, value) in f_eq
                        #println(key, " ", value)
                        bool_add += (getNC_var(fin, key,false).==value)
                        bCounter+=1
                    end
                    # Look for >
                    for (key, value) in f_gt
                        bool_add += (getNC_var(fin, key,false).>value)
                        bCounter+=1
                    end
                    # Look for <
                    for (key, value) in f_lt
                        bool_add += (getNC_var(fin, key,false).<value)
                        bCounter+=1
                    end
    
                    # If all were true, bool_add woule be bCounter!
                    idx = findall(bool_add.==bCounter)
                    ProgressMeter.next!(p; showvalues = [(:File, a), (:N_pixels, size(idx))])
                    # Read data only for non-empty indices
                    if length(idx) > 0
                        #print(size(lat_in_))
                        mat_in =  zeros(Float32,(length(lat_in_[:,1]),length(dGrid)))
                        dim = size(mat_in)
                        # Read in all entries defined in JSON file:
                        co = 1
                        for (key, value) in dGrid
                            #println(key,", ",value)
                            mat_in[:,co]=getNC_var(fin, value,false)
                            co += 1
                        end
    
                        iLat_ = ((lat_in_[idx,:].-latMin)/(latMax-latMin)*length(lat)).+1
                        iLon_ = ((lon_in_[idx,:].-lonMin)/(lonMax-lonMin)*length(lon)).+1
    
                        favg_all!(mat_data, mat_data_variance, mat_data_weights, compSTD, iLat_,iLon_,mat_in[idx,:],length(idx),dim[2],nGrid, latMin, latMax, lonMin,lonMax, length(lat), length(lon), points )
                        #println("Read ", a, " ", length(idx))
                    else
                        #println("Read ", a, " ", length(idx))
                    end
                    close(fin)
                else
                    ProgressMeter.next!(p; showvalues = [(:File, a), (:N_pixels, 0)])
                    close(fin)
                end

                     
            catch
               println("Error in file caught")
            end
        end # End of loop through all files
            
        # Filter all data, set averages, 
        dims = size(mat_data)
        println("Averaging final product...")
        if maximum(mat_data_weights)>0
            #println("mat_data_weights:", size(mat_data_weights))
            out_mat_n[cT,:,:] = mat_data_weights
            co = 1
            for (key, value) in dGrid
                da = round.(mat_data[:,:,co],sigdigits=6)
                da[mat_data_weights.<1e-10].= FillValue
                out_mat_av[cT,:,:,co] = da #
                if compSTD
                    da = round.(sqrt.(mat_data_variance[:,:,co] ./ mat_data_weights)  ,sigdigits=6)
                    da[mat_data_weights.<1e-10].= FillValue
                    out_mat_std[cT,:,:,co] = da 
                end
                co += 1
            end
        else

    
        end
        fill!(mat_data,0.0)
        fill!(mat_data_weights,0.0)
        fill!(mat_data_variance,0.0)
        end # end of loop over eachindex(dates)

    # Define optional fields conditionally
    optional_fields = if compSTD
        (std_data = out_mat_std,) # Note the trailing comma for a single-element Named Tuple
    else
        # Return an empty Named Tuple if the condition is false
        # This allows the merge (splatting) in the final return 
        (;) 
    end

    println("Done.")
    
    return (
        lat        = lat,
        lon        = lon,
        time       = dates,
        averaged_data = out_mat_av,
        n = out_mat_n,
        # ... and any other fields ...
        optional_fields... # Splatting merges the key-value pairs from the optional Named Tuple
    )
end # run_gridding function


function run_gridding_km(files_any::Any;
        json_dict::String = "jsonFiles/S5_O3_dict.json",
        target_km::Float64 = 10.0, # The desired resolution in km
        monthly::Bool=false, 
        compSTD::Bool=false, 
        latMin::Float64=-90.0, 
        latMax::Float64=90.0, 
        lonMin::Float64=-180.0, 
        lonMax::Float64=180.0,
        startDate_str::String="2025-11-25", 
        stopDate_str::String="2015-11-26",
        dDays::Int64=8, 
        oversample_temporal::Float64=1.0,
        FillValue=-999
    )
    # Note: difference to original version when calling this function is that dLat/dLon are not needed anymore as both are now controlled by target_km
    
    println("Running S5 gridding (Python mode, in-memory, sinusoidal km-scale ragged grid)")

    # Convert PyList{Any} or Any other list to a Julia Vector{String}
    # at the same time convert list into a Set for fast lookup 
    files_set = Set(convert(Vector{String}, collect(files_any)))
    
    startDate = DateTime(startDate_str)
    stopDate = DateTime(stopDate_str)
    dDay  = monthly ? Dates.Month(dDays) : Dates.Day(dDays)
    println(startDate, " ", stopDate)
    
    dates = collect(startDate:dDay:stopDate)
    nTime = length(dates)

    # Define Latitude Grid (Constant distance in degrees)
    # degree latitude is approximately 111.32 km
    dLat = target_km / 111.32
    lat_centers = collect((latMin + dLat/2):dLat:(latMax - dLat/2))
    nLat = length(lat_centers)

    # Pre-calculate Longitude Bins per Row to maintain target_km
    # n_lons_per_row[i] = how many bins in the i-th latitude row
    lon_span = lonMax - lonMin
    n_lons_per_row = Int[]
    for l in lat_centers
        # Earth circ in km at this lat
        circ_at_lat = 40075.0 * cos(deg2rad(l))
        # 2. km distance of span at this latitude (km)
        regional_dist_km = circ_at_lat * (lon_span / 360.0)
        # 3. Number of bins needed to cover that distance at target_km resolution
        # use max(1, ...) to ensure at least one bin exists
        push!(n_lons_per_row, max(1, ceil(Int, regional_dist_km / target_km)))
    end

    #**************************************
    # 3. Parse JSON and setup
    jsonDict = JSON.parsefile(json_dict)
    d2 = jsonDict["basic"]
    dGrid = jsonDict["grid"]
    nVars = length(dGrid)
    
    # hard-coded value to divide each polygon (sub-pixel splitting), can be changed:
    nGrid = 10;
    #println("nGrid is set to: ", nGrid)
    #####
    points = zeros(Float32, (nGrid, nGrid, 2))

    # Read all filters:
    f_eq = getFilter("filter_eq",jsonDict)
    f_gt = getFilter("filter_gt",jsonDict)
    f_lt = getFilter("filter_lt",jsonDict)

    # Get file naming pattern (needs YYYY MM and DD in there)
    fPattern = jsonDict["filePattern"]
    # Get main folder for files:
    # folder   = jsonDict["folder"] # not needed as files are given via python... files already exist

    
    #**************************************

    # Time counter
    cT = 1
    p1 = Progress(cT) 
    
    # Initialize output containers
    # We use a Vector (time) of Vectors (lat) of matrices (lon x variables)
    # This replaces the old 4D Array: out_mat[time, lon, lat, var]
    final_output_av  = [] # Structure: final_output_av[time][lat_row][lon_bin, var_idx]
    final_output_n   = [] # Structure: final_output_n[time][lat_row][lon_bin]
    final_output_std = compSTD ? [] : nothing # Only if needed. Structure: final_output_std[time][lat_row][lon_bin, var_idx]

    
    # Process Time Steps
    for cT in eachindex(dates)
        d = dates[cT]
        
        println("\n Gridding time slice $d ($cT/$nTime)")
        ProgressMeter.next!(p1; showvalues = [(:Time, d)])
    
        # Initialize temporary output arrays (time step specific):
        mat_data     = [zeros(Float32, (n_lons_per_row[i], nVars)) for i in 1:nLat]
        mat_weights  = [zeros(Float32, n_lons_per_row[i]) for i in 1:nLat]
        # mat_n        = [zeros(Float32, (n_lons_per_row[i])) for i in 1:nLat] # not needed, can use mat_weights directly
        # Initialize variance if needed
        mat_variance = compSTD ? [zeros(Float32, (n_lons_per_row[i], nVars)) for i in 1:nLat] : nothing
        
        # Define the search window
        # This matches the range from orgiginal code: [d, d + dDay*oversample_temporal - 1 Day]
        date_range = d:Dates.Day(1):d + dDay * oversample_temporal - Dates.Day(1)

        ### find files within this time step from provided list
        # Create specific patterns for each date
        target_patterns = [replace(fPattern, 
                            "YYYY" => Dates.format(di, "yyyy"), 
                            "MM"   => Dates.format(di, "mm"), 
                            "DD"   => Dates.format(di, "dd")) 
                           for di in date_range]
        
        # Create a single Regex that combines all dates: e.g., (pattern1|pattern2|pattern3)
        # replace '*' with '.*' to make it a valid regex wildcard
        combined_regex = Regex(join(replace.(target_patterns, "*" => ".*"), "|"))
        
        # Filter once
        files_tmp = filter(f -> occursin(combined_regex, f), files_set)
        # Convert the Set to a Vector (Array) so it supports positional and boolean indexing.
        files_vector = collect(files_tmp)        
        
        println("  -> Files found (", length(files_vector), "): ")#, files)
        #######

        # check if file is not empty
        fileSize = Int[];
        for f in files_vector # Iterate over the Vector
            fileSize = [fileSize;stat(f).size]
        end

        println("File size: ", fileSize) 
        
        # Loop through all files within this time step
        p = Progress(length(files_vector))

        for a in files_vector[fileSize.>0]
            try
                fin = Dataset(a)
                println("Reading: $a")
 
                lat_center = getNC_var(fin, d2["lat"], false)
                lon_center = fixlon(getNC_var(fin, d2["lon"], false))

                # --- convert 0–360 to –180–180 (center coords) if needed ---
                lon_center = fixlon(lon_center)
                
                # Quick check if file overlaps region
                bool_found = (lat_center.>latMin) .+ (lat_center.<latMax) .+ (lon_center.>lonMin) .+ (lon_center.<lonMax)
                # If a sounding passes all criteria, its value is 4
                
                if any(bool_found.==4)
                    
                    lat_in_ = getNC_var(fin, d2["lat_bnd"], true)
                    lon_in_ = fixlon(getNC_var(fin, d2["lon_bnd"], true))

                    # Transpose if orders are swapped
                    dim = size(lat_in_)
                    if dim[1]==4
                        lat_in_ = lat_in_'
                        lon_in_ = lon_in_'
                    end
    
                    # --- convert 0–360 to –180–180 (bounds) if needed ---
                    lon_in_ = fixlon(lon_in_)

                    # Find all indices within lat/lon bounds:
                    minLat = minimum(lat_in_, dims=2)
                    maxLat = maximum(lat_in_, dims=2)
                    minLon = minimum(lon_in_, dims=2)
                    maxLon = maximum(lon_in_, dims=2)
 
                    # Get indices within the lat/lon bounding box and check filter criteria (the last one filters out data crossing the date boundary):
                    bool_add = (minLat[:,1].>latMin) .+ (maxLat[:,1].<latMax) .+ (minLon[:,1].>lonMin) .+ (maxLon[:,1].<lonMax) .+ ((maxLon[:,1].-minLon[:,1]).<50)
                    # If a sounding passes all 5, its value in bool_add is 5.
                    
                    ## Filtering 
                    bCounter = 5
                    # Look for equalities
                    for (key, value) in f_eq
                        #println(key, " ", value)
                        bool_add += (getNC_var(fin, key,false).==value)
                        bCounter+=1
                    end
                    # Look for >
                    for (key, value) in f_gt
                        bool_add += (getNC_var(fin, key,false).>value)
                        bCounter+=1
                    end
                    # Look for <
                    for (key, value) in f_lt
                        bool_add += (getNC_var(fin, key,false).<value)
                        bCounter+=1
                    end
                    # ...
                    
                    idx = findall(bool_add .== bCounter)
                    ProgressMeter.next!(p; showvalues = [(:File, a), (:N_pixels, size(idx))])

                    if !isempty(idx)
                        # Load data variables
                        mat_in = zeros(Float32, (length(lat_in_[:,1]), nVars))
                        for (co, (key, value)) in enumerate(dGrid)
                            mat_in[:, co] = getNC_var(fin, value, false)
                        end

                        # Call the KM-scale splatting function
                        favg_all_km!(mat_data, mat_variance, mat_weights, compSTD, 
                                     lat_in_[idx, :], lon_in_[idx, :], mat_in[idx, :], 
                                     length(idx), nVars, nGrid, latMin, latMax, 
                                     lonMin, lonMax, nLat, n_lons_per_row, points)
                    else
                    ProgressMeter.next!(p; showvalues = [(:File, a), (:N_pixels, 0)])
                    close(fin)
                    end    
                end
                close(fin)
                next!(p)
            catch e
                println("Error in file $a: $e")
            end
        end # loop over files in this time step
        
        # ************************
        # Postprocess this time step
        # Handle FillValues and store result for this time step
        # Since it's a ragged array, process row by row
        for r in 1:nLat
            for c in 1:n_lons_per_row[r]
                if mat_weights[r][c] < 1e-10
                    mat_data[r][c, :] .= FillValue
                end
            end
        end
        # 2. Post-processing variance
        if compSTD
            for r in 1:nLat
                for c in 1:n_lons_per_row[r]
                    if mat_weights[r][c] >= 1e-10
                        for z in 1:nVars
                            # Convert M2 (sum of squares) to standard deviation
                            # σ = sqrt( M2 / weight )
                            mat_variance[r][c, z] = round(sqrt(mat_variance[r][c, z] / mat_weights[r][c]), sigdigits=6)
                        end
                    else
                        mat_variance[r][c, :] .= FillValue
                    end
                end
            end
        end
        
        # Push to a output list
        push!(final_output_av, mat_data)
        push!(final_output_n, mat_weights)
        if compSTD
            push!(final_output_std, mat_variance)
        end
    end # end of loop over eachindex(dates)

    
    # Define optional fields conditionally
    optional_fields = if compSTD
            (std_data = final_output_std,)
        else
            (;) 
        end
    
    println("Gridding Complete.")

    return (
    lat_centers = lat_centers,
    n_lons_per_row = n_lons_per_row,
    time = dates,
    averaged_data = final_output_av,
    n = final_output_n,    
    optional_fields...
    )

end # run_gridding_km function
 


end # module