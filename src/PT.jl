# This file contains the parts of PRZMTools.jl that are used by cranberry
"""
Converts a single "NxRate@Interval" string into two arrays of rates and
intervals. It splits on the "x" and the @ characters. If there is a "kg"
in the Rate expression, it assumes the application rate is given in kg,
otherwise it assums it's in grams and converts to kg.

Returns a 2-tuple of 1-D arrays, one for the application rate and the
second for the interval from the previous application.

Set includeinterval to *true* if this is to be appended to another schedule.
this starts counting the intervals from *interval* rather than zero """
function decode1app(s, includeinterval = false)
  x = split(s,'x'; limit=2)
  naps, y = if size(x,1) == 1 (1, x[1]) else (parse(Int, x[1]), x[2]) end
  z = rsplit(y,'@'; limit=2)
  interval = if size(z,1) > 1  # there was an interval after the "@"
    parse(Int, replace(strip(z[2]), r"[^0-9].*"=>""))
  elseif (naps == 1 | !includeinterval) # interval not needed, set to -1
    -1
  else
    @warn "Need interval, assuming 7 days"
    7
  end
  unitmult = occursin("kg", z[1]) ? 1.0 : 0.001
  rate = parse(Float64, replace(strip(z[1]), r"[^0-9.].*"=> "")) * unitmult
  rates = repeat([rate], inner=[naps])
  offsetstart = includeinterval ? interval : 0
  offset = range(offsetstart, step=interval, length=naps) |> collect
  return(of = offset, rt=rates)
end

function decodeapp(s::String)
  appstrs = split(s, ",")
  apparrays = [decode1app(appstrs[i], i>1) for i in 1:size(appstrs,1)] #need to set includeinterval true after first call
  for i in 2:size(apparrays,1)
    apparrays[i].of .+= apparrays[i-1].of[end]
  end
  return(of=vcat(map(x->x.of, apparrays)...), rt = vcat(map(x->x.rt, apparrays)...))
end


# Water body
@with_kw struct WaterBody{R<:Real} @deftype R # lines are for vvwmtransfer/PWC files
  burial::Bool = false # line 34/54
  PrBen = -0.5 # Distribution of pesticide on eroded sediment. use negative to choose "varying"
  Ddx = 1e-8        # Ln 39/59, mass xfer coefficient, in m/day
  bendep = 0.05     # Ln 41/61 (in m)
  porosity = 0.5    # Ln 42/62
  BD = 1.35         # Ln 43/63
  FrOCW = 0.04      # Ln 50/70
  FrOCB = 0.04      # Ln 44/64
  DOCW = 5.0        # Ln 51/71
  DOCB = 5.0        # Ln 45/65
  BNMAS = 0.006     # Ln 46/66, Benthic biomass
  BfAC = 1.19       # Ln 47/67, photolysis parameter DFAC
  SUSED = 30.0      # Ln 48/68, suspended solids
  CHL = 0.005       # Ln 49/69, Chlorophyll; affects photolysis attenuation only
  PLMAS = 0.4       # Ln 52/72, water column biomass
  Area = 10000.0    # Ln 36/56, 1 ha in m^2
  Depth = 0.8       # Ln 37/57, initial water depth, in m
  MaxDepth = 0.8    # Ln 38/58, max water depth, in m
  FieldA = 100000.0 # Ln 35/55
  FlowAvgDays::Int = 30  # Days for flow averaging; use negative number for no flow averaging
  baseflow = 0.0 # not saved in .PWC file
  # Ω = 1e-8   # = Ddx/bendep
end

"""
    PWCWB(lines::Array{String,1})

Read a WaterBody from a PWC split into an array of lines
"""
function PWCWB(lines::Array{String,1})
#  lines = readlines(pwcfile)
  burial = parse(Bool, lowercase(lines[54]))
  prbvary, prbvalue = split(lines[60],',')[1:2]
  PrBen = lowercase(prbvary) == "true" ? -parse(Float64,prbvalue) : parse(Float64, prbvalue)
  floavg,dummy,flowdays = split(lines[53],',')
  if flowdays =="" flowdays = "30" end # put in a value if there is none, should be the case only when floavg==false
  FlowAvgDays = lowercase(floavg) == "true" ? parse(Int,flowdays) : -parse(Int,flowdays)
  # return DataFrame(lineNo = [59,61,62,63,70,64,71,65,66,67,
  #                                            68,69,72,55,56,57,58,51],
  #                   value = lines[[59,61,62,63,70,64,71,65,66,67,68,69,72,55,56,57,58,51]] )
  values = [parse(Float64,x) for x in lines[[59,61,62,63,70,64,71,65,66,67,
                                             68,69,72,56,57,58,55]]]
  WaterBody(burial, PrBen, values..., FlowAvgDays,  0.0) #  baseflow not saved in PWC file, but in GUI
end

"""
    PFAMWB(lines::Array{String,1})

Convert lines 55-90 of a PFAM input file to a water body file. Uses the
larges *fill* value in the PFAM file to set the depth"""
function PFAMWB(lines::Array{String,1})
  vals = parse.(Float64, lines[vcat(21,23:30,35,36)])
  fillstr = split(lines[1],',')
  fill = maximum(parse.(Float64,fillstr[1:findfirst(fillstr.=="")-1]))
  WaterBody(false, -0.5, vals[1:6]..., vals[9], 0.0, 0.0, # DOCB and BNMAS not in PFAM, set 0
            vals[10], vals[7], vals[8], 0.0,  # PLMAS not in PFAM
            10000.0, fill, fill, 1000000.0, 30, 0.0)
end
#wb = WaterBody("../testing/TestChemTest.PWC")

"""
    WaterBody(fn::String)

Read water body data from either a PWC or PFAM input file, or return a standard
water body. If fn is "IR", "P", or "W", then return the Index Reservoir, EPA farm 
pond or PMRA wetland, respectively. Otherwise assumes fn is a file name for a 
PFAM or PWC file, and attempts to read it. Water body parameters that are not in
the file are assigned default values"""
function WaterBody(fn)
  if fn == "IR" # set up the EPA Index Reservoir
    WaterBody(FieldA = 172.8e4, Area = 5.26e4, Depth=2.74, MaxDepth=2.74)
  elseif fn == "P" # EPA pond
    WaterBody( Depth=2.0, MaxDepth=2.0, FlowAvgDays=0.0)
  elseif fn == "W" # PMRA 80cm wetland
    WaterBody(FlowAvgDays = 0.0)
  else # assume fn is a filename
    lines = readlines(fn)
    if occursin("PFAM", lines[1])
      PFAMWB(lines[55:90])
    elseif occursin("PWC", lines[1])
      PWCWB(lines) # PWCWB was originally written to read the file
    else
      @error "The file $fn must be either a PWC or a PFAM input file"
    end
  end
end



# Pesticide
@with_kw struct Pesticide @deftype Vector{Float64}
  Name::String = ""
  Code::String = ""
  Sf::Matrix{Float64} = [-0.099 0.0 ; 0.059 -0.058]
  nchm::Int64 = size(Sf,1)
  SQ10::Float64 = 2.0
  STref::Float64 = 25.0
  Wf::Matrix{Float64} = zeros(nchm, nchm)
  WQ10::Float64 =  2.0
  WTref::Float64 = 25.0
  Bf::Matrix{Float64} = zeros(nchm, nchm)
  BQ10::Float64 =  2.0
  BTref::Float64 = 25.0
  Pf::Matrix{Float64} = zeros(nchm, nchm)
  Platref::Float64 = 35.0
  Hf::Matrix{Float64}= zeros(nchm, nchm)
  HQ10::Float64 = 2.0
  Kd = zeros(nchm)
  Koc = zeros(nchm)
  Kclay = zeros(nchm)
  frN = ones(Float64, nchm) # Freundlich 1/n
  VP = zeros(nchm) .+ 1e-10
  VPtemp = zeros(nchm) .+ 25.0
  Sol = zeros(nchm) .+ 100000.0
  DifAir = zeros(nchm) .+ 5000.0
  HofH = zeros(nchm) .+ 50000.0
  MWT = zeros(nchm) .+ 100.0
  PlantUptake = zeros(nchm)
  PlantVol = zeros(nchm)
  PlantDecay::Matrix{Float64} = zeros(nchm, nchm)
  PlantFolEx = zeros(nchm) .+ 0.5
  @assert size(Sf) == size(Wf) == size(Bf) == size(Pf) == size(Hf) == size(PlantDecay)
  @assert size(Sf,1)  == nchm
end

"""Split on commas up to a number of values (nchm) and convert to
Vector{Float64}. Converts blank values to zero, and should always return
a vector of length nchm"""
function parsePWCline(line::String, nchm::Int)
  x = split(line, ',', limit=nchm+1)
  x = x[1:min(size(x,1),nchm)]  # drop any element more than nchm
  [y == "" ? 0.0 : parse(Float64, y) for y in x]
end
parsePWCline(F::IOStream, nchm::Int) =  parsePWCline(readline(F))

"Read a PWC .SWI or .PWC file into a Pesticide"
function PWCChem(lines::Array{String,1}; vptemp = 25.0, freundlichN=1.0)
  chemname = lines[2]
  nchm = parse(Int,lines[3])
  iskoc = parse(Bool, lines[4] |> lowercase)
  lines = [parsePWCline(lines[i],nchm) for i in 5:28]
  #close(F)
  wf = pzconv(lines[2]..., lines[15][1:nchm-1]...)
  bf = pzconv(lines[4]..., lines[16][1:nchm-1]...)
  sf = pzconv(lines[9]..., lines[19][1:nchm-1]...)
  pf = pzconv(lines[6]..., lines[17][1:nchm-1]...)
  hf = pzconv(lines[8]..., lines[18][1:nchm-1]...)
  ff = pzconv(lines[11]..., lines[20][1:nchm-1]...)
  q10 = lines[end][1]
  zro = zeros(Float64, nchm)
  kd, koc = iskoc ? (zro, lines[1]) : (lines[1], zro)
  vpT = zro .+ vptemp # not in PWC input
  frN = zro .+ freundlichN
  Pesticide(chemname,"",sf,nchm,q10,lines[10][1],wf,q10,lines[3][1], bf,
            q10, lines[5][1],pf,lines[7][1], hf, q10, kd, koc, zro, frN,
            lines[13], vpT, lines[14],lines[21], lines[23], lines[12],
            zro, zro, zeros(Float64,nchm,nchm), zro .+ 0.5)
end

function PFAMChem(lines::Array{String,1},name, vptemp=25.0, freundlichN=1.0)
  nchm = parse(Int, lines[5])
  plines = [parsePWCline(lines[i],nchm) for i in 6:25]

  wf = pzconv(plines[1]..., plines[16][1:nchm-1]...)
  bf = pzconv(plines[2]..., plines[17][1:nchm-1]...)
  sf = pzconv(plines[3]..., plines[18][1:nchm-1]...)
  pf = pzconv(plines[4]..., plines[19][1:nchm-1]...)
  hf = pzconv(plines[5]..., plines[20][1:nchm-1]...)
  Pesticide(Name = name, nchm=nchm, Sf = sf, Wf=wf, Bf=bf,Pf=pf, Hf=hf,
            MWT = plines[6], VP=plines[7],Sol = plines[8], Koc=plines[9],
            WTref=plines[10][1], BTref=plines[11][1],STref=plines[12][1],
            Platref=plines[13][1], HofH = plines[14])
end

"""Read Pesticide data from either a PWC or PFAM input file. Any
parameters not in the file are assigned default values"""
function Pesticide(fn; vptemp=25.0, freundlichN=1.0)
  lines = readlines(fn)
  if occursin("PFAM", lines[1])
    chemname = fn[1:end-3] # not in PFAM file, use filename
    PFAMChem(lines, chemname, vptemp, freundlichN)
  elseif occursin("PWC", lines[1])
    PWCChem(lines; vptemp=vptemp, freundlichN=freundlichN) # PWCWB was originally written to read the file
  else
    @error "The file $fn must be either a PWC or a PFAM input file"
  end
end

"convert rate to halflife without changing zeros"
ratehalflife(r) = (r == 0.0 ? 0.0 : log(2.0)/r)

"convert a single halflife or rate constant to a 1x1 matrix with a rate"
function pzconv(P::Float64, halflife=true) # one float
  zeros(Float64,1,1) .- (halflife ? ratehalflife(P) : P)
end

"""Convert parent halflife, daugher halflife and conversion fraction
to a 2x2 fractional transfer coefficient matrix"""
function pzconv(P::Float64, D::Float64, FF::Float64, halflife = true) # three floats
  if halflife
    p = ratehalflife(P)
    return([-p 0.0; p*FF -ratehalflife(D)])
  else
    return [-P 0.0; P*FF -D]
  end
end

function pzconv(P::Float64, D::Float64, G::Float64, ffPD::Float64,
                ffDG::Float64, halflife=true) # five floats
  if halflife
    p,d,g = ratehalflife.((P, D, G))
    return [-p 0 0 ; p*ffPD -d 0.0 ; 0.0 d*ffDG -g]
  else
    return [-P 0 0 ; P*ffPD -D 0 ; 0 D*ffDG -G]
  end
end

"""Returns two vectors: one of halflives (or rates) and the other of formation
fractions. for a 3x3 matrix, a parent→granddaughter formation factor is
returned (middle element of second returned vector). This is not used by PWC,
but can be used by PRZM"""
function pzconv(f::Matrix{Float64}, halflife=true)
  nchm = size(f,1)
  if nchm==1
    # make sure returned rate is positive
    return ([abs(halflife & (f[1,1] != 0.0) ? log(2)/f[1,1] : f[1,1])], [0.0])
  elseif nchm == 2
    ff = f[1,1] == 0.0 ? 0.0 : -f[2,1]/f[1,1]
    if halflife
      p,d = -ratehalflife.(diag(f))
      return([p,d],[ff])
    else
      return (-diag(f), [ff])
    end
  elseif nchm==3
    ffPD = f[1,1] == 0.0 ? 0.0 : -f[2,1]/f[1,1]
    ffDG = f[2,2] == 0.0 ? 0.0 : -f[3,2]/f[2,2]
    ffPG = f[3,3] == 0.0 ? 0.0 : -f[3,1]/f[1,1]
    if halflife
      p, d, g = -ratehalflife.(diag(f))
      ([p,d,g], [ffPD, ffPG, ffDG])
    else
      return (-diag(f), [ffPD, ffPG, ffDG])
    end
  else
    @warn "pzconv given a matrix not of size 1,2 or 3 using as if 3x3"
    x = f[1:3, 1:3]
    ffPD = x[1,1] == 0.0 ? 0.0 : -x[2,1]/x[1,1]
    ffDG = x[2,2] == 0.0 ? 0.0 : -x[3,2]/x[2,2]
    ffPG = x[3,3] == 0.0 ? 0.0 : -x[3,1]/x[1,1]
    if halflife
      p, d, g = -ratehalflife.(diag(x))
      return([p,d,g], [ffPD, ffPG, ffDG])
    else
      return (-diag(x), [ffPD, ffPG, ffDG])
    end
  end
end