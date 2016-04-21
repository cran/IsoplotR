focus <- function(){
    cur.plot <- grDevices::dev.cur()
    if (cur.plot>1){
        grDevices::dev.off()
        grDevices::dev.new()
    }
}
