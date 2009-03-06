setClass("NTimes",
  representation(
    n="integer",
    cases="integer",
    bt="integer",
    et="integer"
  )
)

nTimes <- function(ntimes) {
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)
  out <- new("NTimes",
    n=as.integer(ntimes),
    cases=as.integer(lt),
    bt=as.integer(bt),
    et=as.integer(et)
  )
  out
}
