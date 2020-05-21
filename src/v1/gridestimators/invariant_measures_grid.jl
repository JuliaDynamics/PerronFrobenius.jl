
function invariantmeasure(pts, method::SingleGrid)
    tog = transopergenerator(pts, method)
    to = tog()
    invariantmeasure(to)
end
