var membersFile = new XMLHttpRequest();
membersFile.open("GET", "version.txt", false);
membersFile.send();
lines = membersFile.responseText.split("\n")

// assuming it's a github.io page the first part of the FQDN is the username on
// github and the first part of the path (after the initial /) is the repo name
user = window.location.host.split('.')[0]
repo_name = window.location.pathname.split('/')[1]
text='<img src="https://github.com/'+user+'/'+repo_name+'/actions/workflows/CI.yml/badge.svg" style="display:block;margin-left: auto;margin-right: auto;">';

lines.slice().sort((a, b) => a - b).reverse().forEach(element => {
    if (element != '') { // skips empty element at end of file
        text+='<a href="index_'+element+'.html">Build #'+element+'</a>'
    }
});
    
    

sidebar=document.getElementsByClassName('sidebar');
sidebar[0].innerHTML=text;
