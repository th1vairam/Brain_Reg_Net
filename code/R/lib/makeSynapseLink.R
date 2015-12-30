makeSynapseLink <- function(fileId, parentId, linkName = NULL, annotations = list()){
  library(synapseClient)
  
  # Get synapse properties of the source object
  obj = synapseClient::synGet(fileId, downloadFile=F)
  if (is.null(linkName))
    linkName = obj@properties$name
  
  # Create a linksTo list
  linksTo = list(
    targetVersionNumber = obj@properties$versionNumber,
    targetId = obj@properties$id
  )
  
  # Create a linkClass list
  linkClass = list(
    concreteType  = 'org.sagebionetworks.repo.model.Link',
    linksTo	= linksTo,
    entityType = 'org.sagebionetworks.repo.model.Link',
    parentId	= parentId,
    name = linkName,
    annotations	= annotations,
    linksToClassName = 'org.sagebionetworks.repo.model.FileEntity')
  
  synapseClient::synRestPOST('/entity',linkClass)
}