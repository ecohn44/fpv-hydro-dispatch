using MAT
using Printf
using CSV
using DataFrames
using Dates
using Base.Filesystem
using LaTeXStrings
include("functions.jl")

function fullsim_dataload()

    LMP_path = string("data/fullsim/price-data.csv");
    inflow_path = string("data/fullsim/lake-powell-inflow-daily.csv");
    solarrad_path = "data/fullsim/solar-rad-capacity.csv";
    release_path = string("data/fullsim/lake-mead-release.csv");
    mead_storage_path = string("data/fullsim/lake-mead-storage-daily.csv");
    powell_storage_path = string("data/fullsim/lake-powell-storage-daily.csv");

    # Read in Local Marginal Price for Lake Mead
    RTP = DataFrame(CSV.File(LMP_path)); # currently just 2022
    # datetime handling
    split_datetime = split.(RTP.OPR_DT, "/")
    # Create new columns for Month, Day, and Year by extracting each part
    RTP.Month = getindex.(split_datetime, 1)  # Extract Month
    RTP.Day   = getindex.(split_datetime, 2)  # Extract Day
    RTP.Year  = getindex.(split_datetime, 3)  # Extract Year

    # Read in system inflow (lake powell inflow)
    inflow = DataFrame(CSV.File(inflow_path));
    # convert to m3/s from cfs
    inflow.inflowm3 = cfs_to_m3s(inflow.inflow);
    # Filter for 2022-2023
    inflow_s = subset(inflow, :datetime => dt -> endswith.(dt, "22") .|| endswith.(dt, "23"))
    # datetime handling
    split_datetime = split.(inflow_s.datetime, "/")  # Apply split to each row
    # Create new columns for Month, Day, and Year by extracting each part
    inflow_s.Month = getindex.(split_datetime, 1)  # Extract Month
    inflow_s.Day   = getindex.(split_datetime, 2)  # Extract Day
    inflow_s.Year  = getindex.(split_datetime, 3)  # Extract Year

    # Read in mead storage levels 
    mstorage = DataFrame(CSV.File(mead_storage_path));
    # convert from ac-ft to m3
    mstorage.storagem3 = af_to_m3(mstorage.storage);
    # Filter for 2022-2023
    mstorage_s = subset(mstorage, :datetime => dt -> endswith.(dt, "22") .|| endswith.(dt, "23"))
    mstorage_s.Month = getindex.(split_datetime, 1)  # Extract Month
    mstorage_s.Day   = getindex.(split_datetime, 2)  # Extract Day
    mstorage_s.Year  = getindex.(split_datetime, 3)  # Extract Year

    # Read in powell storage levels 
    pstorage = DataFrame(CSV.File(powell_storage_path));
    # convert from ac-ft to m3
    pstorage.storagem3 = af_to_m3(pstorage.storage);
    # Filter for 2022-2023
    pstorage_s = subset(pstorage, :datetime => dt -> endswith.(dt, "22") .|| endswith.(dt, "23"))
    pstorage_s.Month = getindex.(split_datetime, 1)  # Extract Month
    pstorage_s.Day   = getindex.(split_datetime, 2)  # Extract Day
    pstorage_s.Year  = getindex.(split_datetime, 3)  # Extract Year

    # Read in Water Contract Data
    release = DataFrame(CSV.File(release_path));
    # convert from ac-ft to m3
    release.releasem3 = af_to_m3(release.Result);
    # Filter for 2022-2023
    release_s = subset(release, :Datetime => dt -> endswith.(dt, "22") .|| endswith.(dt, "23"))
    release_s.Month = getindex.(split_datetime, 1)  # Extract Month
    release_s.Day   = getindex.(split_datetime, 2)  # Extract Day
    release_s.Year  = getindex.(split_datetime, 3)  # Extract Year

    # Read in solar radiation capacity
    alpha = DataFrame(CSV.File(solarrad_path));
    for col in names(alpha)
        alpha[!, col] = min_max_normalize(alpha[!, col])  # Normalize and replace the data
    end
    select!(alpha, Not(:Column1))

    # Combine Daily data
    daily = DataFrame(
        year = inflow_s.Year,
        month = inflow_s.Month,
        day = inflow_s.Day,
        inflow = inflow_s.inflowm3,
        release = release_s.releasem3,
        storage = pstorage_s.storagem3 + mstorage_s.storagem3
    )

    return daily, alpha, RTP
end