model_reactions <- function(rxn_id,rxn_name,reactions,genes,subsystem){
  reactions <- sapply(reactions, rxn2sybil)
  rever <- as.vector(sapply(reactions, reversibility))
  comp <- as.vector(sapply(reactions, react_compartment))
  ub<-rep(1000,length(rever))
  lb<-rep(-1000,length(rever))
  lb[rever=="irreversible"]<-0
  model_reac<-as.data.frame(cbind(as.character(rxn_id),as.character(rxn_name),reactions,rever,as.vector(comp),lb,ub,0,as.character(genes),as.character(subsystem)),stringsAsFactors=FALSE)
  colnames(model_reac)<-c("abbreviation","name","equation","reversible","compartment","lowbnd","uppbnd","obj_coef","rule","subsystem")
  return(model_reac)
  }
