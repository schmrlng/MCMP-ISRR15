import Base.show
show(io::IO, x::NearNeighborCache) = print(io, "$(length(x)) samples")
show(io::IO, x::PointRobot2D) = print(io, "Point Robot 2D")
show(io::IO, x::MPSolution) = print(io, "$(x.status): $(x.cost)")
show(io::IO, x::RealVectorMetricSpace) = print(io, "$(x.dist) Metric Space")
show(io::IO, x::LQGPath) = print(io, "Path of length $(length(x.path))")
