# MCMP-ISRR15

Code associated with [Monte Carlo Motion Planning for Robot Trajectory Optimization Under Uncertainty](http://arxiv.org/abs/1504.08053) submitted to [ISRR 2015](http://www.isrr-2015.org/). This code relies on [MotionPlanning.jl](https://github.com/schmrlng/MotionPlanning.jl). This [notebook](http://nbviewer.ipython.org/github/schmrlng/MCMP-ISRR15/blob/master/MCMP%20Example.ipynb) provides some examples of visualizing MCMP algorithm progress.

## Installation Instructions
These installation instructions assume that you're starting from a blank slate system and are looking to run [```MCMP Example.ipynb```](http://nbviewer.ipython.org/github/schmrlng/MCMP-ISRR15/blob/master/MCMP%20Example.ipynb) locally on your computer. If, for example, you have a python installation that you'd prefer to use instead of Anaconda, or otherwise have some pieces of the tech stack already installed, you may have to do some error-message chasing to add missing dependencies.
- Install [Anaconda](https://store.continuum.io/cshop/anaconda/)
- Install [Julia](http://julialang.org/downloads/)
- Open a Julia terminal window and run the following commands

    ```julia
    Pkg.clone("https://github.com/schmrlng/MotionPlanning.jl.git")
    Pkg.add("IJulia")
    Pkg.checkout("Interact")      # master branch of Interact is required for working with IPython 3
    Pkg.build("IJulia")
    ```

- Clone this repository
- Launch IJulia, a browser-based notebook interface for Julia, by either typing

    ```julia
    using IJulia
    notebook()
    ```

  at the ```julia>``` prompt or running ```ipython notebook --profile julia```from the command line.
- Navigate to the folder containing this repository in the IJulia browser window and open ```MCMP Example.ipynb```.
- File an [issue](https://github.com/schmrlng/MCMP-ISRR15/issues/) if the above instructions don't work and I'll try to help you out!
