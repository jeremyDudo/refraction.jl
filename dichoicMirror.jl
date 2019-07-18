using LinearAlgebra
using Plots, ProgressMeter

# Fresnel Coefficients "D" matrix either in TE or TM polarization
function De(n, theta)
    return [ 1 1 ; n*cos(theta) -n*cos(theta)]
end
function Dm(n, theta)
    return [cos(theta) -cos(theta) ; n n ]
end

function kFunc(n, λ, theta)
    """
    n: between ~1-10
    lambda: in nanometers (400-700)
    theta: angle for particular medium (in degrees)
    """
    # c0 = 2.998e+17;                     # nm/s
    return (n*2*pi)/(λ) * cos(theta)
end

# Propagation matrix
function P(n,d,λ,theta)
    k = kFunc(n,λ,theta)
    return [ exp(k*d*im) 0 ; 0 exp(-k*d*im)]
end


function Mirror(n1, d1, n2, d2, nSub, dSub, ntot, polarization)
    theta1 = asin(sin(45)/n1);
    theta2 = asin(n1*sin(theta1)/n2);

    if polarization == "TE"
        D = De
    elseif polarization == "TM"
        D = Dm
    else 
        println("Polarization must be either \"TE\" or \"TM\"")
    end

    # be very careful with number of layers here, can break
    # case 1 should have an even number of layers
    # case 2 should have an odd number of layers
    if nSub == n1
        return P(nSub, dSub) * ( inv(D(n1,d1)) * D(n2,d2) * P(n2, d2) * inv(D(n2,d2)) * D(n1,d1) * P(n1, d1) )^((ntot - 1)/2) * D(1,45)
    
    elseif nSub == n2
        return P(nSub, dSub) * inv(D(n2,d2)) * D(n1,d1) * P(n1, d1) * ( inv(D(n1,d1)) * D(n2,d2) * P(n2, d2) * inv(D(n2,d2)) * D(n1,d1) * P(n1, d1) )^((ntot - 2)/2) * D(1,45)
    else
        println("nSub must equal n1 or n2")
        return 0
    end
end

Mirror(50, 0.005, 90, 0.005, 90, 20, 9, "TE")

struct MirrorParams
    # usually between 1-10
    n1::Float16
    n2::Float16

    # in nanometers
    d1::Float64
    d2::Float64

    # larger than nanometers
    dSub::Float16

    # number of layers
    totalLayers::Int

    # TE or TM
    Polarization::String
end

function mirr(mp::MirrorParams, λ)
    theta1 = asin(sin(45)/mp.n1);
    theta2 = asin(mp.n1*sin(theta1)/mp.n2);

    if mp.Polarization == "TE"
        D = De
    elseif mp.Polarization == "TM"
        D = Dm
    else 
        println("Polarization must be either \"TE\" or \"TM\"")
    end

    # be very careful with number of layers here, can break
    # case 1 should have an even number of layers
    # case 2 should have an odd number of layers
    if mp.totalLayers//2 == 0
        nSub = mp.n1
        return inv(D(nSub, mp.dSub))*D(1, 45) * P(nSub, mp.dSub) * ( inv(D(mp.n1,mp.d1)) * D(mp.n2,mp.d2) * P(mp.n2, mp.d2) * inv(D(mp.n2,mp.d2)) * D(mp.n1,mp.d1) * P(mp.n1, mp.d1) )^((mp.totalLayers)/2) * inv(D(1,45))*D(mp.n1,theta1)
    else
        nSub = mp.n2
        return inv(D(nSub, mp.dSub))*D(1, 45) * P(nSub, mp.dSub) * inv(D(mp.n2,mp.d2)) * D(mp.n1,mp.d1) * P(mp.n1, mp.d1) * ( inv(D(mp.n1,mp.d1)) * D(mp.n2,mp.d2) * P(mp.n2, mp.d2) * inv(D(mp.n2,mp.d2)) * D(mp.n1,mp.d1) * P(mp.n1, mp.d1) )^((mp.totalLayers - 1)/2) * inv(D(1,45))*D(mp.n1,theta1)
    end
end

function mirr2(n1,n2,d1,d2,dSub,totalLayers,Polarization,λ)
    theta1 = asin(sin(45)/n1);
    theta2 = asin(n1*sin(theta1)/n2);
    # println(theta1, theta2)

    if Polarization == "TE"
        D = De
    elseif Polarization == "TM"
        D = Dm
    else 
        println("Polarization must be either \"TE\" or \"TM\"")
    end

    # be very careful with number of layers here, can break
    # case 1 should have an even number of layers
    # case 2 should have an odd number of layers
    if totalLayers//2 == 0
        nSub = n1
        thetaSub = theta1         
        return inv(D(nSub, dSub))*D(1, 45) * P(nSub, dSub, λ, thetaSub) * ( inv(D(n1,d1)) * D(n2,d2) * P(n2, d2, λ, theta2) * inv(D(n2,d2)) * D(n1,d1) * P(n1, d1, λ, theta1) )^((totalLayers)/2) * inv(D(1,45))*D(n1,theta1)
    else
        nSub = n2
        thetaSub = theta2 # inv(D(nSub, dSub))*D(1, 45) *
        return inv(D(nSub, dSub))*D(1, 45) * P(nSub, dSub, λ, thetaSub) * inv(D(n2,d2)) * D(n1,d1) * P(n1, d1, λ, theta1) * ( inv(D(n1,d1)) * D(n2,d2) * P(n2, d2, λ, theta2) * inv(D(n2,d2)) * D(n1,d1) * P(n1, d1, λ, theta1) )^((totalLayers - 1)/2) * inv(D(1,45))*D(n1,theta1)
    end
end

function reflect2(n1,n2,d1,d2,dSub,totalLayers,Polarization,λ)
    matr = mirr2(n1,n2,d1,d2,dSub,totalLayers,Polarization,λ)
    #  println(matr)
    return ((norm(matr[3]/matr[4])))^2
end

x = 400:0.1:700;

# works:
# 1.46, 4.6, 52, 52, 1000, 21, "TE", x

# not so much 
# 1.46, 2.2, 95, 95, 100, 100, "TM", x
y = reflect2.(1.46, 4.6, 52, 52, 1000, 21, "TM", x);
y[194]  # == 100% 
y[73]   # == 0%

plot(x,y, ylims=(-0.1,1.1))


