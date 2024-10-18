var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ECCO","category":"page"},{"location":"#ECCO","page":"Home","title":"ECCO","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ECCO.jl Julia package.","category":"page"},{"location":"","page":"Home","title":"Home","text":"ECCO stands for Estimating the Circulation and Climate of the Ocean. \nECCO is a 40+ years running effort initiated by Pr. Carl Wunsch at MIT, and which now involves distributed academic projects accross the US (MIT, NASA JPL, UCSD, UTA, WHOI), Europe (AWI, UH), and the UK (BAS).","category":"page"},{"location":"","page":"Home","title":"Home","text":"this package is in early development stage","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ECCO, ECCO.toy_problems, ECCO.glacier_model, ECCO.Lorenz_models]","category":"page"},{"location":"#ECCO.toy_problems.enzyme_ex1-Tuple{}","page":"Home","title":"ECCO.toy_problems.enzyme_ex1","text":"ECCO.toy_problems.enzyme_ex1()\n\nusing ECCO\n(f,f_ad,x)=ECCO.toy_problems.enzyme_ex1()\nf(x)\nf_ad(x)\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.toy_problems.enzyme_ex2-Tuple{}","page":"Home","title":"ECCO.toy_problems.enzyme_ex2","text":"ECCO.toy_problems.enzyme_ex2()\n\nusing ECCO\n(f,f_ad,x)=ECCO.toy_problems.enzyme_ex2()\nf(x...)\nf_ad(x...)\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.toy_problems.enzyme_ex3-Tuple{}","page":"Home","title":"ECCO.toy_problems.enzyme_ex3","text":"ECCO.toy_problems.enzyme_ex3()\n\nusing ECCO\n(f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex3()\nf(x,y)\nf_ad(x,y)\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.toy_problems.enzyme_ex4-Tuple{}","page":"Home","title":"ECCO.toy_problems.enzyme_ex4","text":"ECCO.toy_problems.enzyme_ex4()\n\nusing ECCO\n(f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex4()\nf(x,y)\nf_ad(x,y)\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.toy_problems.optim_ex1-Tuple{}","page":"Home","title":"ECCO.toy_problems.optim_ex1","text":"ECCO.toy_problems.optim_ex1()\n\nusing ECCO\n(f,x0,x1,result)=ECCO.toy_problems.optim_ex1()\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.toy_problems.optim_ex2-Tuple{}","page":"Home","title":"ECCO.toy_problems.optim_ex2","text":"ECCO.toy_problems.optim_ex2()\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.toy_problems.optim_ex3-Tuple{}","page":"Home","title":"ECCO.toy_problems.optim_ex3","text":"ECCO.toy_problems.optim_ex3()\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.glacier_model.forward_problem-Tuple{Array, Int64, Vararg{Float64, 4}}","page":"Home","title":"ECCO.glacier_model.forward_problem","text":"forward_problem(xx::Array, nx::Int, dx::Float64, xend::Float64, \n\t\t\t\tdt::Float64, tend::Float64)\n\nSimple, 1D mountain glacier model inspired from the book Fundamentals of Glacier Dynamics,  by CJ van der Veen, and which was translated to Julia by S Gaikwad.\n\nSee https://sicopolis.readthedocs.io/en/latest/AD/tutorial_tapenade.html#mountain-glacier-model\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.glacier_model.integrate-Tuple{}","page":"Home","title":"ECCO.glacier_model.integrate","text":"integrate()\n\nV=ECCO.glacier_model.integrate()\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.Lorenz_models.L63-Tuple{}","page":"Home","title":"ECCO.Lorenz_models.L63","text":"L63(; nt=10000)\n\nSee https://en.wikipedia.org/wiki/Lorenz_system\n\nusing ECCO, CairoMakie\nx,y,z=ECCO.Lorenz_models.L63()\nlines(x,y)\n\n\n\n\n\n","category":"method"},{"location":"#ECCO.Lorenz_models.L96-Tuple{}","page":"Home","title":"ECCO.Lorenz_models.L96","text":"L96(; N=5, F=8)\n\nSee https://en.wikipedia.org/wiki/Lorenz96model\n\nusing ECCO, CairoMakie\nstore=ECCO.Lorenz_models.L96()\nlines(store[1,:]); lines!(store[2,:]); lines!(store[end,:])\ncurrent_figure()\n\n\n\n\n\n","category":"method"}]
}
