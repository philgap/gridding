module S5Grid

using Base, Dates, Printf
using NCDatasets
using Statistics
using Glob
using JSON
using ProgressMeter

export run_gridding

### some helper functions...
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

# Still need to make sure the corners are read in properly!
function getNC_var(fin, path, DD::Bool)
    loc = split(path ,r"/")
    #@show loc
    if length(loc)==1
        return fin[path].var[:]
    elseif length(loc)>1
        gr = []
        for i in 1:length(loc)-1
            if i==1
                gr = fin.group[loc[i]]
            else
                gr = gr.group[loc[i]]
            end
        end
	   #@show gr
        #println(loc[end])
        si = size(gr[loc[end]])
        #dimnames = gr[loc[end]].dimnames
        #println(dimnames)

        # DD means there is a 2nd index for footprint bounds of dimension 4!
        if DD
            if si[1]==4
                return reshape(gr[loc[end]].var[:],4,prod(si[2:end]))'
            elseif si[end]==4
                return reshape(gr[loc[end]].var[:],prod(si[1:end-1]),4)
            end
        else
	#@show gr[loc[end]]
            return reshape(gr[loc[end]].var[:],prod(si))
        end
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

# Convert 0–360° to -180–180° if needed
fixlon(x) = ifelse.(x .> 180, x .- 360, x)


function run_gridding(files_any::Any;
        json_dict::String = "jsonFiles/S5_O3_dict.json",
        monthly::Bool=false, 
        compSTD::Bool=false, # only false possible for now...
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

    println("Running S5 gridding (Python mode, in-memory)")
    
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
    folder   = jsonDict["folder"] # not needed as files are given via python... files already exist
    
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
    out_mat_n  = zeros(Float32,(length(dates),length(lon),length(lat),length(dGrid)))
    out_mat_std = compSTD ? zeros(Float32,(length(dates),length(lon),length(lat),length(dGrid))) : nothing # This is needed only if compSTD is true

    for cT in eachindex(dates)
        d = dates[cT]
    
        println("Gridding time slice ", d, " (",cT,"/", nTime,")")
        ProgressMeter.next!(p1; showvalues = [(:Time, d)])
    
        # 1. Define the search window
        # This matches the range from orgiginal code: [d, d + dDay*oversample_temporal - 1 Day]
        date_range = d:Dates.Day(1):d + dDay * oversample_temporal - Dates.Day(1)
    
        # 2. Generate all target filenames (full paths) within the window
        target_files = String[]
        for di in date_range
            
            # Step 2a: Substitute the date components into the file pattern
            # The filePattern now generates the exact filename/glob string we are looking for
            filePattern = replace(fPattern, 
                "YYYY" => @sprintf("%04i", Dates.year(di)),
                "MM"   => @sprintf("%02i", Dates.month(di)),
                "DD"   => @sprintf("%02i", Dates.day(di))
            )
    
            # Step 2b: need a way to combine the folder and the pattern for matching
            # Since we are comparing against a list of full paths in all_files,
            # generate the exact pattern that would match the full path.
            # This is where the old `glob` behavior is slightly tricky to replicate.
            
            # need the full file path pattern, e.g., "/data/nc_files/*_G_V_20250101*.nc"
            full_path_pattern = joinpath(folder, filePattern)
            
            # collect all possible patterns (e.g., one for 20250101, one for 20250102, etc.)
            push!(target_files, full_path_pattern)
        end
    
        # 3. Filter the original list (files_set) based on the target patterns
        # use a custom filter function here because the target_files contain '*' wildcards.
        # check if each existing file matches any of the target patterns.
    
        files = filter(files_set) do existing_file_path
            # Check if the existing file matches any pattern in target_files
            return any(target_pattern -> occursin(Regex(replace(target_pattern, "*" => ".*")), existing_file_path), target_files)
        end

        # Convert the Set to a Vector (Array) so it supports positional and boolean indexing.
        files_vector = collect(files)
        
        # Example output:
        println("  -> Files found (", length(files_vector), "): ")#, files)

        fileSize = Int[];
        for f in files_vector # Iterate over the Vector
            fileSize = [fileSize;stat(f).size]
        end

        println("File size: ", fileSize) 
        
        # Loop through all files
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

end # module
