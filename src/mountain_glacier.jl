
module glacier_model

using Enzyme

@inline function update_h(h::Float64, b::Float64)
	if h < b
		h = b
	end
	return h
end

"""
    forward_problem(M0=0.004; dt=1/12.0, nt=6*5000, dx = 1.0, nx = 30)

Simple, 1D mountain glacier model inspired from the book Fundamentals of Glacier Dynamics, 
by CJ van der Veen, and which was translated to Julia by S Gaikwad.

See https://sicopolis.readthedocs.io/en/latest/AD/tutorial_tapenade.html#mountain-glacier-model

```
V=glacier_model.forward_problem(0.002)
```
"""
function forward_problem(M0=0.004; dt=1/12.0, nt=6*5000, dx = 1.0, nx = 30)
	rho = 920.0
	g = 9.2
	n = 3
	A = 1.e-16
	C = 2*A/(n+2)*(rho*g)^n*(1.e3)^n
	bx = -0.0001
	#M0 = .004
	M1 = 0.0002
	tend = nt*dt

	D = zeros(nx)
	phi = zeros(nx)

	xx = zeros(nx+1)
	xarr = [(i-1)*dx for i in 1:(nx+1)]

	M = zeros(nx+1)
	b = zeros(nx+1)
	M .= M0 .- xarr .* M1 .+ xx
	b .= 1.0 .+ bx .* xarr

	h = zeros(nx+1)
	h_capital = zeros(nx+1)

	h .= ones(nx+1) .* b
	h_capital .= h .- b

	for t in 1:nt
		D .= C .* ((h_capital[1:nx] .+ h_capital[2:nx+1]) ./ 2.0).^(n+2) .* ((h[2:nx+1] .- h[1:nx]) ./ dx).^(n-1)
		phi .= -D .* (h[2:nx+1] .- h[1:nx]) ./ dx

		h[2:nx] .= h[2:nx] .+ M[2:nx] .* dt .- dt/dx .* (phi[2:nx] .- phi[1:nx-1])
#		h[2:nx] .= update_h.(h[2:nx], b[2:nx])
		h[2:nx] .= [max(h[i],b[i]) for i in 2:nx]

		h_capital .= h .- b
	end

	V = sum(Array(h_capital[1:nx+1].*dx))
	return V
end

adjoint_problem(x)=autodiff(Reverse, forward_problem, Active(x[1]))

end

