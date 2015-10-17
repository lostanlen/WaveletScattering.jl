subscripts = (PathKey(:time), PathKey(:chunk))
x = rand(Float32, 1024, 16)
node = FourierNode(rand(Float32, 1024, 10), 1, subscripts)
