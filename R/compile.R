#' @export
compile.acre <- function(dev = FALSE, debug.mode = FALSE){
  if(dev){
    dir = paste0(getwd(), '/inst/TMB/')
  } else {
    dir = paste0(system.file(package = "acre"), "/TMB/")
  }
  
  if(file.exists(paste0(dir, "acreTmb.o"))) unlink(paste0(dir, "acreTmb.o"))
  if(file.exists(paste0(dir, "acreTmb.dll"))) unlink(paste0(dir, 'acreTmb.dll'))

  if(!debug.mode){
    TMB::compile(paste0(dir, "acreTmb.cpp"), framework = "CppAD")
  } else {
    if(Sys.info()['sysname'] == 'Windows'){
      TMB::compile(paste0(dir, "acreTmb.cpp"), "-O1 -g", DLLFLAGS="")
    } else {
      TMB::compile(paste0(dir, "acreTmb.cpp"), "-O0 -g")
    }
    
  }

}
