.First.lib<-function(lib,pkg){
    require("rama",  quietly=TRUE) || stop("rama package not found")
    library.dynam("bridge",pkg,lib)
}
