  $('referencesTarget').append( $('.references') );
  $( init );
  
  function init() {
  
   // Move the paragraph from #myDiv1 to #myDiv2 and #myDiv3
    $('.referencesTarget').append( $('.references>p') );
  }
