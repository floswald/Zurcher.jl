# Zurcher.jl

* This is a teaching-oriented package for the Bus Engine Replacement Model after Rust (1987) used for my [Computational Econ Course](https://floswald.github.io/NumericalMethods/)
* We implement both _naive_ nfxp and mpec estimation. It is _naive_ nfxp because it uses standard VFI instead of the much more performant polyalgorithm developed by Rust in his paper. For a thorough benchmarking exercise I refer to 
    1. ECTA Comment by [Fedor Iskhakov, Jinhyuk Lee, John Rust, Bertel Schjerning, Kyoungwon Seo](https://www.econometricsociety.org/publications/econometrica/2016/01/01/comment-“constrained-optimization-approaches-estimation)
    1. Matlab implementation which includes the polyalgorithm and analytic derivatives for likelihood function distributed as part of the [DSE2019](https://github.com/dseconf/DSE2019/tree/master/02_DDC_SchjerningIskhakov/code/zurcher) summer school. Several parts of my code have been copied and modified from that code base.

### What is the point of this package?

* The main point is to demonstrate the relative easiness with which we can tackle an MPEC problem with the [`JuMP.jl`](https://jump.dev) package.
* JuMP is like AMPL, but for free and embedded in a proper programming language. 
* For a similar, even more impressive demonstration of this please visit [https://github.com/UBCECON567/BLPDemand.jl](https://github.com/UBCECON567/BLPDemand.jl)



```@autodocs
Modules = [Zurcher]
```


end
