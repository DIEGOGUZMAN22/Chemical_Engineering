function RK4(ODEfun,tspan,y0,h)
    n = (xf - x0)/h + 1
    y = zeros(Int(n), length(y0))
    for j = 1:length(y0)
        y[1,j] = y0[j]
    end
    x = collect(x0:h:xf)
    for j = 1:length(x)-1
        k1 = zeros(neq)
        k2 = zeros(neq)
        k3 = zeros(neq)
        k4 = zeros(neq)
        for i = 1:neq
            k1[i] = f(x[j],y[j,:])[i]
            k2[i] = f(x[j] + (1/2) * h, y[j,:] + (1/2) * k1[:] * h)[i]
            k3[i] = f(x[j] + (1/2) * h, y[j,:] + (1/2) * k2[:] * h)[i]
            k4[i] = f(x[j] + h, y[j,:] + k3[:] * h)[i]
            y[j+1,i] = y[j,i] + (1/6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * h
        end

    end
    return x,y
end
