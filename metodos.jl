unction euler(f, t₀, y₀, tₙ; n = 100)
    h = (tₙ - t₀) / n
    #t = linspace(t₀, tₙ, n + 1)
    t = range(t₀, stop=tₙ, length=n+1)
    y = zeros(n + 1)
    y[1] = y₀
    for i = 1:n
        y[i+1] = y[i] + f(t[i], y[i]) * h
    end

    return t, y
end

function euler_aperfeicoado(f, t₀, y₀, t; n = 100)
    h = (tₙ - t₀)/n
    #t = linspace(t₀, tₙ, n + 1)
    t = range(t₀, stop=tₙ, length=n+1)
    y = zeros(n + 1)
    y[1] = y₀
    α = β = 0.5
    γ = δ = 1
    for i = 1:n
        k₁ = f(t[i],y[i])
        k₂ = f(t[i] + δ * h, y[i] + k1 * γ * h)
        y[i+1] = y[i] + (α * k1 + β * k2) * h
    end

    return t, y
end

function midpoint(f, t₀, y₀, t; n = 100)
    h = (tₙ - t₀)/n
    #t = linspace(t₀, tₙ, n + 1)
    t = range(t₀, stop=tₙ, length=n+1)
    y = zeros(n + 1)
    γ = δ = 0.5
    α = 0
    β = 1
    y[1] = y₀
    for i=1:n
        k1 = f(t[i], y[i])
        k2 = f(t[i] + h * δ,y[i] + γ * h * k1)
        y[i+1] = y[i] + (α * k1 + β * k2) * h
    end

    return t,y
end

function rungekutta4(f, t₀, y₀, tₙ; n=100)
    h = (tₙ - t₀) / n
    y[1] = y₀
    y = zeros(n + 1)
    #t = linspace(t₀, tₙ, n + 1)
    t = range(t₀, stop = tₙ, length = n + 1)
    for i=1:n
        k1 = f(t[i], y[i])
        k2 = f(t[i] + h * 0.5, y[i] + h * k * 0.5)
        k3 = f(t[i] + h * 0.5, y[i] + h * k2 * 0.5)
        k4 = f(t[i] + h, y[i]+ k3*h)
        y[i+1] = y[i] + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)
    end
    return t, y
end

function heun(f, t₀, y₀, tₙ; n = 100)
  h = (tₙ - t₀) / n   
  y[1] = y₀
  #t = linspace(t₀, tₙ, n + 1)
  t = range(t₀, stop=tₙ, length=n+1)
  y = zeros(n + 1)
  δ = γ = 2/3
  α = 0.25
  β = 0.75
  for k = 1:n
    k₁ = f(t[k], y[k])
    k₂ = f(t[k] + δ * h, y[k] + γ * h * k₁)
    y[k + 1] = y[k] + h * (α * k₁ + β * k₂)
  end
  return t, y
end
