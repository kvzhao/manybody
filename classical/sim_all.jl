include("./ising_model.jl")
using HDF5
#using PyPlot

function thermalization(s::Array{Int8}, L::Int64, T::Float64, THERMO)
    for i in 1:THERMO
        update(s, L, T)
    end
end

# Though we dont like global
const J = 1.0
const TMAX = 3.5
const TMIN = 1.5
const NUM_SIM = 101
const NUM_BINS = 10
# Total simulation time = MCSteps x Bin size
# N = BxM
const THERMO = 20000
const MCSTEPS = 10000
const fMCS = Float64(MCSTEPS)

## Simulation Informations
for L in [8,16,32,64,128]
    N = L^2
    # Output 
    outfile = "L$L.h5"
    println(outfile)
    ## Simulation over different temperature
    TList = linspace(TMAX, TMIN, NUM_SIM)
    spins = init_state(L)
    # Results
    absM_T=[]
    M_T=[]
    E_T=[]
    M2_T=[]
    E2_T=[]
    M4_T=[]
    E4_T=[]

    println("Start 2D Ising Monte Carlo Simulation")
    tic()
    for T in TList
        thermalization(spins, L, T, THERMO)

        # all: N 
        absM_all=zeros(NUM_BINS)
        M_all=zeros(NUM_BINS)
        M2_all=zeros(NUM_BINS)
        M4_all=zeros(NUM_BINS)
        E_all=zeros(NUM_BINS)
        E2_all=zeros(NUM_BINS)
        E4_all=zeros(NUM_BINS)

        for bin in 1:NUM_BINS
            absMbin = 0.0
            Mbin = 0.0
            M2bin = 0.0
            M4bin = 0.0
            Ebin = 0.0
            E2bin = 0.0
            E4bin = 0.0
            for step in 1:MCSTEPS
                update(spins, L, T)
                # measurement
                # TODO: optimize this
                E = cal_eng(spins, L)
                M = cal_mag(spins)
                # energy and magnetization
                eng = E/N
                mag = M/N
                absMbin += abs(mag)
                Mbin += mag
                M2bin += mag^2
                M4bin += mag^4
                Ebin += eng
                E2bin += eng^2
                E4bin += eng^4
            end
            absM_all[bin] = absMbin/fMCS
            M_all[bin] = Mbin/fMCS
            M2_all[bin] = M2bin/fMCS
            M4_all[bin] = M4bin/fMCS
            E_all[bin] = Ebin/fMCS
            E2_all[bin] = E2bin/fMCS
            E4_all[bin] = E4bin/fMCS
        end

        # Calcualte mean and deviation
        absM_mean = mean(absM_all)
        M_mean = mean(M_all)
        M2_mean = mean(M2_all)
        M4_mean = mean(M4_all)
        E_mean = mean(E_all)
        E2_mean = mean(E2_all)
        E4_mean = mean(E4_all)

        append!(absM_T, absM_mean)
        append!(M_T, M_mean)
        append!(M2_T, M2_mean)
        append!(M4_T, M4_mean)
        append!(E_T, E_mean)
        append!(E2_T, E2_mean)
        append!(E4_T, E4_mean)

        println("$T: E: $E_mean, M: $M_mean, |M|: $absM_mean")#, Cv: $Cv, X: $X")
    end
    toc()
    println("Done.")

    # Save the result
    h5open(outfile, "w") do file
        write(file, "T", Array{Float64}(TList))
        write(file, "absM", Array{Float64}(absM_T))
        write(file, "M", Array{Float64}(M_T))
        write(file, "M2", Array{Float64}(M2_T))
        write(file, "M4", Array{Float64}(M4_T))
        write(file, "E", Array{Float64}(E_T))
        write(file, "E2", Array{Float64}(E2_T))
        write(file, "E4", Array{Float64}(E4_T))
    end
    println("Save the data to $outfile")
end