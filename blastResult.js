
function finalSubmit(dbIsNt, targetForm, inputName, submitFormName, dbType){
  var formArr=targetForm.split(":");
  var idArray=new Array();
  if(dbType == 0){
      if(dbIsNt==1){
          document.forms[submitFormName].db.value="Nucleotide";
      }
      else{
          document.forms[submitFormName].db.value="Protein";
      }
  }
  for(j=0; j<formArr.length; j++){
    for(var i=0; i<document.forms[formArr[j]].elements.length; i++){
      var theElem=document.forms[formArr[j]].elements[i];
      if(theElem.type=="checkbox"&&theElem.name==inputName&&theElem.checked){
        if(!isIdIn(theElem.value, idArray)){
     
         idArray[idArray.length]=theElem.value;
         }
      }
    }
    if(dbType == 0){
        document.forms[submitFormName].term.value="";
        for(i=0; i<idArray.length; i++){  
            if(i ==0 ) {
                document.forms[submitFormName].term.value += idArray[i];
            } else {
                document.forms[submitFormName].term.value += "," + idArray[i];
            }
	}	
    } else if (dbType == 1){
        document.forms[submitFormName].val.value="";
        for(i=0; i<idArray.length; i++){    
            if(i ==0 ) {
                document.forms[submitFormName].val.value += idArray[i];
            } else {   
                document.forms[submitFormName].val.value += "," + idArray[i];
            }
	}	
    }
  }
    document.forms[submitFormName].submit();
}

function isIdIn(id, idArray){
  var idSeen=false;

  for(i=0; i<idArray.length; i++){
    if(id==idArray[i]){
      idSeen=true;
      break;
    }
  }
  return idSeen;
}

function handleCheckAll(mode, targetForm, inputName){
 
  var formArr=targetForm.split(":");
  for(j=0; j<formArr.length; j++){
    for(var i=0; i<document.forms[formArr[j]].elements.length; i++){
      var theElem=document.forms[formArr[j]].elements[i];
      if(theElem.type=="checkbox"&&theElem.name==inputName){
        if(mode=="select"){
          theElem.checked=true;
        }
        else if(mode=="deselect"){
          theElem.checked=false;
        }
      }
    }
  }
 
}

function synchronizeCheck(id, formName, inputName, isChecked){
 
  for(var i=0; i<document.forms[formName].elements.length; i++){
    var theElem=document.forms[formName].elements[i];
    if(theElem.type=="checkbox"&&theElem.name==inputName&&id==theElem.value){    
        theElem.checked=isChecked; 
       
    }
  }

}


function uncheckable(formName, inputName){
  for(var i=0; i<document.forms[formName].elements.length; i++){
    var theElem=document.forms[formName].elements[i];
    if(theElem.type=="checkbox"&&theElem.name==inputName){    
        theElem.checked=0; 
       
    }
  }

}
