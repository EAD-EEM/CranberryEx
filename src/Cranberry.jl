module Cranberry
using Parameters, DataFrames, Dates, CSV
using LinearAlgebra: eigen, diagm, diag
include("PT.jl") 
export cranSch, WBMat, soildeg, cranberryyear, PWCsoil, PWCChem, WaterBody

#Q10TimeAdjust(T; T0 = 20.0, Q10=2.0) = Q10.^((T0 .- T)./10.0) |>  cumsum

# """Generates an array with columns of application rate and days after
# first flood, given an initial application date, an application schedule
# string and the day of the flood"""
# function appsched(firstapp::Date, appsched::String, firstflood::Date)
#   appoffset = Dates.value(firstapp - firstflood)
#   amt, offset = PRZMTools.decodeapp(appsched)
#   offset = offset .+ appoffset
#   hcat(amt, offset)
# end

"""
  cranSch(weafile::String, app1date::Date, appsch::String, flooddates, T0=20.0, Q10=2.0)
  cranSch(wea::DataFrame, app1date::Date, appsch::String, flooddates, T0=20.0, Q10=2.0)

Creates a cranberry schedule DataFrame containing the time related inputs for the 
cranberry model.

# Arguments
- `weafile` is a file name of a PRZM .dvf/.met or .wea formatted weather file

- `wea` is a dataframe with two columns, for date and temperature

- `app1date is the date of the first application

- `appsch` is a string of the form `nxrate@d` defining the use pattern of the pesticide

- `flooddates` is an array of Dates giving the start day for each flood. The last
  date is the date the final bog is drained

# Details
The retured DataFrame contains columns for day, date, temperature, temperature-adjusted
day, application amount and a boolean where *true* gives flood dates.
"""
function cranSch(weafile::String, app1date::Date, appsch::String, flooddates, T0=20.0, Q10=2.0)
F = open(weafile, "r")
line = readline(F) # for determining file type
seekstart(F)
wea = if occursin(',',line) # this is assumed to be a .wea (comma delimited) file
  w = CSV.File(F, header=false) |> DataFrame
  w.Date = Date.(w.Column3, w.Column1, w.Column2)
  w[!,[9,6]] # keep date and temperature only
else  #assume fixed width (read as multiple space delimited) .dvf or .met file
  w = DataFrame!(CSV.File(F, dateformat="mmddyy", ignorerepeated=true, 
            header=false, delim=" ", types=Dict(1=>Date)))
  w[!,1] = map(x->year(x) < 35 ? x+Year(2000) : x+Year(1900), w[!,1])
  w[!,[1,4]]
end
cranSch(wea, app1date, appsch, flooddates, T0, Q10)
end

function cranSch(wea::DataFrame, app1date::Date, appsch::String, flooddates, T0=20.0, Q10=2.0)
  # make the years of application and flooding match the first year of weather
  app1 = Date(year(wea[1,1]), month(app1date), day(app1date))
  floods = map(x->Date(year(wea[1,1]),month(x), day(x)), flooddates)
  @assert wea[1,1] <= app1 "function cranSch: weather file starts $(wea[1,1]), after given date"
  wea = wea[wea[!,1] .>= app1, :]
  wea.Day = 1:nrow(wea)
  wea.adjDay = Q10.^((wea[!,2] .- T0)./10.0) |>  cumsum #Q10TimeAdjust(wea[!,2], T0=T0, Q10=Q10)
  rename!(wea,[:Date, :Temp, :Day, :adjDay])
  wea.app = zeros(Float64, nrow(wea))
  sch = decodeapp(appsch)
  wea[sch.of .+ 1, :app] .= sch.rt
  wea.flood = zeros(Bool, nrow(wea))
  wea.flood[findall(x->x∈floods, wea.Date)] .= true
  wea
end

"""   soluteholdingcapacity(S::WaterBody, koc::Array{Float64,1})

Calculates solute holding capacity for each Koc in the `koc` array and
then returns the total capacity for the water compartment, for the sediment
compartment (including pore water), and for the pore water alone."""
function soluteholdingcapacity(S::WaterBody, koc::Array{Float64,1})
  wvol = S.Depth * S.Area  # volume of water body (vol1 in VVWM)
  m_wsed = S.SUSED * wvol * 0.001 # mass suspended sediment (m_sed_1 in VVWM)
  kd_wsed = koc * S.FrOCW * 0.001
  kd_wbio = 0.436 .* (koc ./ 0.35) .^ 0.907 .* 0.001
  m_wbio = S.PLMAS * wvol * 0.001
  kd_wdoc = koc ./ 0.35 .* 0.000074 # incl factor of 1/1000
  m_wdoc = S.DOCW * wvol * 0.001
  watercap = kd_wsed .* m_wsed .+ kd_wbio .* m_wbio .+ kd_wdoc * m_wdoc .+ wvol
  pvol = S.bendep * S.Area * S.porosity # pore water volume (v2 in VVWM)
  kd_sed = koc .* S.FrOCB * 0.001
  m_sed = S.BD * S.bendep * S.Area * 1000.0
  kd_bio = 0.436 .* (koc ./ 0.35) .^ 0.907 .* 0.001 # VVWM: only use of XKPB and KOW

  m_bio = S.BNMAS * S.Area * 0.001 #sediment, water equiv is m_wbio
  kd_pdoc = koc .* 0.001  # DOC in pore water
  m_pdoc = S.DOCB * pvol * 0.001
  porecap = kd_pdoc .* m_pdoc .+ pvol
  sedcap = kd_sed .* m_sed .+ kd_bio .* m_bio .+ porecap
  (water = watercap, sed = sedcap, pore = porecap)
end

"""    fconv(f)

Convert a matrix which has the rates of outflow from each compartment on
the diagonal into a proper fractional transfer coefficient matrix by replacing
the diagonal with -colSums."""
function fconv(f)
  f1 = deepcopy(f)
  for (i, v) in enumerate(sum(f1, 1))
    f1[i, i] = -v
  end
  f1
end

"""    compartmentEval(f::Array{Float64,2}, strt::Array{Float64,1}, times::Array{Float64,1})

Evaluate a compartment model 

# Arguments
- `f` fractional transfer coefficient matrix *f*, a square matrix 

- `times` The times at which to evaluate the model

- `strt` An array of starting values for each compartment

*f* must be a square matrix, and *f* and strt must have the same dimension.
same as *size(startConc,2)*."""
function compartmentEval(f::Array{Float64,2}, strt::Array{Float64,1}, 
    times::Array{Float64,1})
  # startConc = vcat([1],repeat([0],inner=[size(f1)[1]-1]))
  vals, vecs = eigen(f)
  cnst = vecs \ strt
  c2 = diagm(0 => cnst) * exp.(vals * times')
  result = (vecs * c2)'
end

""" Not used (2021 Feb)
Temperature adjusted evaluation of a compartment model. Adjust temperature
by Q10 method. Assume the temperature vector, *Temp* has one value per time
unit, and returns values at all time units, *i.e.* returns a vector of the
same length as Temp, but is intended to be paired with evenly spaced times
of one time unit."""
function adjEval( f::Array{Float64,2}, strt::Array{Float64,1},
  Temp::Array{Float64,1}, T0, Q10)
  evaltimes = Q10TimeAdujust(Temp, T0, Q10)
  compartmentEval(f, strt, evaltimes)
end

"""Sum the columns of a matrix, and provide a version that just returns
a scalar"""
colsum(x::Array{<:Real,2}) = hcat(sum(x; dims = 1))
colsum(x::Real) = x

#= the next two functions are replacements for f4PDwatsed and WBMat,
   grouping by compartment rather than chemical (i.e. the matrices for
   degradation of a chemical are along the diagonal, and transfers between
   compartments are distributed. Also allows allows an arbitrary
   number of compartments.=#
# String matrices to test makefmatrix
# Pc1 = ["p11-1" "p12-1" "P13-1";"p21-1" "p22-1" "p23-1"; "p31-1" "p32-1" "p33-1"]
# Pc2 = ["p11-2" "p12-2" "P13-2";"p21-2" "p22-2" "p23-2"; "p31-2" "p32-2" "p33-2"]
# Pc3 = ["p11-3" "p12-3" "P13-3";"p21-3" "p22-3" "p23-3"; "p31-3" "p32-3" "p33-3"]
#
# Cp1 = ["x" "c12-1" "c13-1";"c21-1" "y" "c23-1" ; "c31-1" "c32-1" "z"]
# Cp2 = ["x" "c12-2" "c13-2";"c21-2" "y" "c23-2" ; "c31-2" "c32-2" "z"]
# Cp3 = ["x" "c12-3" "c13-3";"c21-3" "y" "c23-3" ; "c31-3" "c32-3" "z"]
#
# ps = [Pc1, Pc2, Pc3]
# cs = [Cp1, Cp2, Cp3]
# np = size(ps,1)
# nc = size(cs,1)
# n = np*nc
# test to see if it's putting stuff where expected:
#ps = [[-3 0; 2 -1], [-2 1;1 -1]]
#cs = [[0 10; 20 0],[0 30; 40 0]]
#makefmatrix(ps,cs)
"""Makes a fractional transfer coefficient matrix for a system with chemical
transformation in multiple compartments. It requires a fractional transfer
coefficient matrix for the chemical transformation in each comparment (a
vector of matrices) and a fractional transfer coefficient matrix for movement
of each chemical between compartments. It assumes there is no loss of chemical
in transfer between compartments and ignores the diagonals on these inputs.

The sizes and numbers of input matrices are
+ deg (list of chemical transfomation matrices): length: number of compartments;
  dimension of matrices: number of chemicals
+ cmp (list of chemical transfers between compartments) : length: number of chemicals;
  dimension of matrices: number of compartments
"""
function makefmatrix(deg, cmp)
  @assert size(unique(size.(deg)), 1) == 1
  @assert size(deg[1], 1) == size(cmp, 1)
  @assert size(cmp[1], 1) == size(deg, 1)
  typ = eltype(eltype(deg))
  if typ <: Real
    out = -hcat(colsum.(deg)...)' #.- vcat(diag.(deg)...)
  end
  np = size(cmp, 1) # number of chemicals, 1 matrix in cmp for each chem
  nc = size(deg, 1) # number of physical compartments, 1 matrix in deg for each
  n = np * nc
  M = zeros(typ, n, n)
  for (i, m) in enumerate(cmp) # copy the transfers between compartments, done first so diagonal can be overwritten when copying the degradation matrices
    M[i:np:end, i:np:end] .= m
  end
  for (i, m) in enumerate(deg)
    k = i * np
    j = k - np + 1
    M[j:k, j:k] .= m
  end
  if typeof(deg[1][1, 1]) <: Real
    for i = 1:size(M, 1)
      M[i, i] = -(out[i] + sum(M[:, i]) - M[i, i])
    end
  end
  M
end

#makefmatrix([Pc1[1:2,1:2],Pc2[1:2,1:2],Pc3[1:2,1:2]], cs[1:2])
#makefmatrix(ps[1:2], [Cp1p2[1:2,1:2],Cp3[1:2,1:2]])

"""    WBMat(bog::WaterBody, pst::Pesticide, temp=-1000.0)
Generate a fraction transfer coefficient matrix from water body and pesticide paraemters

The generated matrix stores the chemical transformations in blocks. With *n* chemicals
the first *n*×*n* block is the sediment compartment and the second is the water 
compartment."""
function WBMat(bog::WaterBody, pst::Pesticide, temp=-1000.0)
  # VVWM has separate k_aer_s and k_anaer_aq, which are set to k_aer_aq and
  # k_anaer_s, respectively, and would multiply (1-fw1) and (1-fw2), below
  if any((pst.Kd .> 0) .& (pst.Koc .<= 0))
    @warn "Cranberry requires a Koc not Kd"
  end
  (Wf, Bf, Hf) = if temp < -500.0  # no temperature adjustment
    (pst.Wf, pst.Bf, pst.Hf)
  else  # adjust water, hydrolysis and sediment degradation matrices
    (pst.Wf .* (pst.WQ10^((temp-pst.WTref)/10.0)),
     pst.Bf .* (pst.BQ10^((temp-pst.BTref)/10.0)), pst.Hf) # hydrolysis no temp adj
  end
  cap = soluteholdingcapacity(bog, pst.Koc)
  # Fraction of solute in the water column
  # orig VVWM: v1=daily_depth*area; v1/(kd_sed_1 * m_sed_1 + kd_bio * m_bio_1 + kd_doc_1 * m_doc_1 + v1)
  fw1 = (bog.Area * bog.Depth) ./ cap.water
  # fraction of solute in the benthic region
  fw2 = (bog.Area * bog.bendep * bog.porosity) ./ cap.sed
  θ = cap.sed ./ cap.water
  γaq = (Wf .+ Hf) .* fw1 .+ Wf .* (1 .- fw1) # water rate;VVWM has k_aq and k_s, but sets k_s=k_aq
  γsed = Bf .* fw2 .+ Bf .* (1 .- fw2)# same for benthic compartment
  Ω = bog.Ddx / bog.bendep * 86400.0 # Ddx/benthic depth, ddx is "mass transfer coefficient" in m/s
  cmp = [[0.0 ωθ; Ω 0.0] for ωθ in Ω .* θ]
  (makefmatrix([γsed, γaq], cmp), cap.pore ./ cap.sed)
end

"""
    soildeg(p::Pesticide, app1::Date, apsch::String, outt::Array{Date,1})
    soildeg(p::Pesticide, sch::DataFrame)

Calculate degradation of a transforming chemical in soil

+ `sch` a DataFrame holding application and temperature data, as returned by *cranSch*
Calculates degradation of a pesticide and its transformation products in soil, returning
an array win one column for each chemical."""
function soildeg(p::Pesticide, app1::Date, apsch::String, outt::Array{Date,1})
  apparray = decodeapp(apsch)  # of = offset from initial date; rt=rate
  outdays = outt .- app1 .|> Dates.value .|> Float64 # day numbers to write results
  result = zeros(Float64, size(outt, 1), p.nchm)
  decaytime = Dates.value(outt[1] - app1) .- apparray.of # time between appl and 1st output
  # solve the system for each application
  conc = zeros(Float64, 1, p.nchm)
  for i in 1:size(apparray.of,1)
    start = zeros(Float64, p.nchm)
    start[1] = apparray.rt[i]
    conc .+= compartmentEval(p.Sf, start, Float64[decaytime[i]])
  end
  if size(outdays,1) > 1
    vcat(conc, compartmentEval(p.Sf, vec(conc), outdays[2:end].-outdays[1]))
  else
    conc
  end
end

function soildeg(p::Pesticide, sch::DataFrame)
#  apparray = decodeapp(sch)  # of = offset from initial date; rt=rate
  floodidx = findall(sch.flood)[1:end-1] # last date ignored becasue it is end only, not transfer
  outdays = sch.adjDay[floodidx]
  # generate a matrix of run times between each app (columns) and each flood
  # (rows). Each row gives the time from an application to each flood. Sum
  # the results over columns to get the total at each flood from all applications.
  decaytimes = floodidx .- permutedims(sch.adjDay[sch.app .> 0.0])
  result = zeros(Float64, size(outdays,1), p.nchm) # values
  apps = sch.app[sch.app .> 0.0]
  for i in size(decaytimes,2)
    start = zeros(Float64, p.nchm)
    start[1] = apps[i]
    #@show apps[i], typeof(r) , r
    result .+= compartmentEval(p.Sf, start, decaytimes[:,i])
  end
  result
end

"""    cranberryharvest(p::Pesticide, wb::WaterBody, sch::DataFrame)

The last flooddate should be the date the final bog is drained, and water
moved to a holding pond."""
function cranberryharvest(p::Pesticide, wb::WaterBody, sch::DataFrame)
  soilconc = soildeg(p, sch)
  floodidx = findall(sch.flood) # index days (sch rows) with flood start;
  result = zeros(Float64,nrow(sch)-findfirst(sch.flood)+1,2*p.nchm)
  residx = floodidx .- (floodidx[1] - 1) # indices into result matching floodidx, but offset
  M,psratio = WBMat(wb, p)
  start = vcat(soilconc[1,:], zeros(Float64, size(soilconc,2)))
  for i in 2:size(floodidx,1)
    # calculate times from the day before the interval
    schrows = floodidx[i-1]:floodidx[i]-1
    adjt = sch.adjDay[schrows] .- sch.adjDay[floodidx[i-1]-1] # adjusted time from start of flood
    result[residx[i-1]:residx[i]-1,:] .= compartmentEval(M, start, adjt)
    if i<size(floodidx,1) # don't set start on last iteration, as there is no soilconc
      start = vcat(soilconc[i,:], result[residx[i]-1,p.nchm+1:end])
    end
  end
  # after last flood, run out to the end of the
  start = result[residx[end]-1,p.nchm+1:end]
  schrows = floodidx[end]:nrow(sch)
  adjt = sch.adjDay[schrows] .- sch.adjDay[floodidx[end]-1] # adjusted time from start of flood
  result[residx[end]:end,p.nchm+1:end] .= compartmentEval(p.Wf, start, adjt)
  result
end

"""    cranberryyear(pst::Pesticide, wb::WaterBody, sch::DataFrame)

Return daily concentrations in sediment and water for *n* cranberry bogs.
Returns a matrix with the sediment concentration for each bog in the first *n*
columns, and one water concentration in the *n+1st* column"""
function cranberryyear(pst::Pesticide, wb::WaterBody, sch::DataFrame)
  @assert sch.app[1] > 0.0  # assumes start on day of first app
  floodidx = findall(sch.flood)
  nfl = size(floodidx,1)-1 # number of floods, last flood date is drain only
  appidx = findall(sch.app .> 0.0)
  @unpack nchm, Sf, Wf = pst
  ɸfrac = 0.5 # fraction of porosity that drains (porosity - field capacity, but FC is not in wb)
  #need nchm * nflood + nchm columns for result
  result = zeros(Float64,size(sch,1),nchm*(nfl+1)) # last set is for water
  start = zeros(Float64, nchm) # amt in water at start of flood cycle
  n = floodidx[end-1] # maybe -1 to be the end of the last flood rather than the first day after it. this is how far we run the dry degradation
  for i in appidx
    start[1] = sch.app[i]
    adjtime = sch.adjDay[i:n] .- (i==1 ?  0.0 : sch.adjDay[i-1])
    result[i:n,1:nchm] .+= compartmentEval(Sf,start,adjtime)
  end
  # copy result to columns for other bogs (they're the same until flooding)
  for i in range(nchm.+1,stop=nchm*nfl,step=nchm)
    result[1:n, i:i+nchm-1] .= result[1:n, 1:nchm]
  end
  # run each flood
  wcol = (nfl)*nchm+1:(nfl+1)*nchm # columns for water concentrations
  drainage = zeros(Float64, nchm)
  for i in 1:size(floodidx,1)-1 # flood events are start-transfer...transfer-end
    scol = (i-1)*nchm+1:i*nchm  # column for soil conc used for iᵗʰ flood
    rcol = vcat(scol, wcol)     # columns used for this flood
    rows = floodidx[i]:floodidx[i+1]-1 # rows for this flood
    # the combined water-sediment evaluation uses non-time-dependent mass-
    # transfer, and time-dependent degradation, so it is run day-by-day:
    psratio = 0.0 # need this after the for loop, so set it up here
    for r in rows
      start = result[r-1, rcol]  # initial concentrations in soil and water
      # insert drainage here, before next flood, each compartment proportional to volume not included because drainage volume will fill pore volume in next bog
      start[nchm+1:end] .+= (drainage .* wb.Depth/(wb.Depth+wb.bendep*ɸfrac)) # flood water
      start[1:nchm] .+= (drainage .* wb.bendep*ɸfrac/(wb.Depth+wb.bendep*ɸfrac)) # into empty pores
      drainage .= 0.0 # add only on flood day
      M, psratio = WBMat(wb, pst, sch[r, :Temp])
      times = Float64[sch.Day[r]-sch.Day[r-1]] # one element array of float for the one time
      result[r, rcol] = compartmentEval(M, start, times) # don't require 1-day time step
    end
    drainage = psratio .* result[floodidx[i+1]-1, scol] # amount in soil pore water
    drainage .*= ɸfrac  # amount in soil water that drains
    # Soil concentration from end of flood to end of run
    adjtime = sch.adjDay[floodidx[i+1]:end] .- sch.adjDay[rows[end]]
    start = result[floodidx[i+1]-1, scol] .- drainage # remove drainage from soil
    result[floodidx[i+1]:end,scol] = compartmentEval(Sf, start, adjtime) # run out soil compartment
  end
  # now run out the water compartment alone (no sediment) simulating a worst-case
  # storage pond (or one with pesticide already in the sediment)
  start = result[floodidx[end]-1, wcol] .+ drainage
  adjtime = sch.adjDay[floodidx[end]:nrow(sch)] .- sch.adjDay[floodidx[end]-1]
  result[floodidx[end]:end,wcol] = compartmentEval(Wf,start,adjtime)
  result
end

end # Module