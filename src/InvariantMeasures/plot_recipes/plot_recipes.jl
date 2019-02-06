import StateSpaceReconstruction:
	minima_and_stepsizes,
    splitaxes,
    Rectangle

using RecipesBase

@recipe function f(r::RectangularInvariantMeasure,
    boxfillfactor::Int = 3, linesegments = true,
    lw = 0.8, lc = :black, lα = 0.5, ls = :do,
    ms = 0.2, mc = :black, mα = 0.3)

    # All bins are assigned a bin origin.
    unique_visited = unique(r.visited_bins_coordinates, dims = 2)

    # Bins may be visited by several times, so we'll get rid
    # of the repetitions.
    n_visited_boxes = size(unique_visited, 2)

    D = size(unique_visited, 1)

    if D != 3 
        throw(ArgumentError("Plot recipe not implemented for dimension $D."))
    end
    

    # Plot the partition grid over the points of the reconstructed
    # orbit.
    legend --> false

    for i = 1:n_visited_boxes
        origin = unique_visited[:, i]
        rect = Rectangle(
            SVector{3, Float64}(unique_visited[:, i]), 
            SVector{3, Float64}(r.ϵ_absolute))

        @series begin 
            rect
        end

        # Get the measure of the box, and fill the box with a
        # number of points scaled to the measure.
        #μ = riv.measure.dist[i]
        #fillpts = fill_hyperrectangle_3D(origin,
        #            stepsizes,
        #            ceil(Int, μ*100)*boxfillfactor)
    end

    @series begin
        seriestype := :scatter3d
        splitaxes(r.points)
    end

end


@recipe function f(r::InducedRectangularInvariantMeasure, plot_boxes = true)

    legend --> false

    # Plot all boxes at the target resolution with nonzero measure.
    bins_with_nonzero_measure_ϵF = r.bins_ϵF[:, findall(r.measure_induced.dist .> 0)]
    n_visited_boxes_ϵF = maximum(size(bins_with_nonzero_measure_ϵF))

    for i = 1:n_visited_boxes_ϵF
        origin = r.bins_ϵF[:, i]
        rect = Rectangle(
            SVector{3, Float64}(bins_with_nonzero_measure_ϵF[:, i]), 
            SVector{3, Float64}(r.ϵF_absolute))

        @series begin 
            lc --> :black
            lw --> 2
            rect
        end
    end

    # Plot all boxes visited at the source resolution.
    n_visited_boxes_ϵj = maximum(size(r.bins_visited_ϵj))

    for i = 1:n_visited_boxes_ϵj
        origin = r.bins_visited_ϵj[:, i]
        rect = Rectangle(
            SVector{3, Float64}(r.bins_visited_ϵj[:, i]), 
            SVector{3, Float64}(r.ϵj_absolute))

        @series begin 
            lc --> :green
            rect
        end
    end

    @series begin
        seriestype := :scatter3d
        mc --> :red
        splitaxes(r.points)
    end
end


@recipe function f(r::AverageRectangularInvariantMeasure, boxes = true, measure = false)

    if boxes && !measure
        #####################
        # TARGET RESOLUTION
        #####################

        # All bins are assigned a bin origin.
        unique_visited = unique(r.measure_ϵF.visited_bins_coordinates, dims = 2)

        # Bins may be visited by several times, so we'll get rid
        # of the repetitions.
        n_visited_boxes = size(unique_visited, 2)
          D = size(unique_visited, 1)

        if D != 3 
            throw(ArgumentError("Plot recipe not implemented for dimension $D."))
        end
        

        # Plot the partition grid over the points of the reconstructed
        # orbit.
        legend --> false

        for i = 1:n_visited_boxes
            origin = unique_visited[:, i]
            rect = Rectangle(
                SVector{3, Float64}(unique_visited[:, i]), 
                SVector{3, Float64}(r.measure_ϵF.ϵ_absolute))

            @series begin 
                lw --> 2
                rect
            end
        end

        ######################################################
        # RESOLUTIONS FROM WHICH THE MEASURE IS INDUCED
        ######################################################

        cols = [:red, :purple, :green, :blue, :yellow, :magenta, :pink]
        for (k, measure) in enumerate(r.induced_measures)
            # Plot all boxes visited at the source resolution.
            n_visited_boxes_ϵj = maximum(size(measure.bins_visited_ϵj))

            for i = 1:n_visited_boxes_ϵj
                origin = measure.bins_visited_ϵj[:, i]
                rect = Rectangle(
                    SVector{3, Float64}(measure.bins_visited_ϵj[:, i]), 
                    SVector{3, Float64}(measure.ϵj_absolute))

                @series begin 
                    lc --> cols[k]
                    rect
                end
            end
        end

        ######################################################
        # THE POINTS FROM WHICH THE MEASURE IS INDUCED
        ######################################################
        @series begin
            seriestype := :scatter3d
            splitaxes(r.measure_ϵF.points)
        end


    end
end





