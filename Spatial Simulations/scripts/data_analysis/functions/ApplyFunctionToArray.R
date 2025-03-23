## ApplyFunctionToArray ## 
## takes a nested list as the genindNestedAray (output of SINS_get_dada_for_analysis, applies a function to each element and returns the output in a list with the same structure as the original nested list ##

ApplyFunctionToArray<-function(nested_list,func,...){
  isNested <- function(l){
    stopifnot(is.list(l))
    for (i in l) {
      if (is.list(i)) return(TRUE)
    }
    return(FALSE)
  }
  
  if(isNested(nested_list)==FALSE){
  
    func_list<-lapply(nested_list,func,...)
  
  }else{
    
   ##create empty nested list###
  func_list<-vector("list",length=length(nested_list))
  
  for (i in 1:length(nested_list)){
    func_list[[i]]<-vector("list",length=length(nested_list[[1]]))
    
    func_list[[i]]<-lapply(nested_list[[i]],func,...)
  }
  }
  
  return(func_list)
}
