

const checkDual=()=>{
    document.getElementById('freq').addEventListener('click',()=>{
        if(parseInt(document.getElementById('freq').value)>3){
            document.getElementsByClassName('ionex_file')[0].style.visibility="hidden"
            document.getElementsByClassName('ionex_file')[0].style.display="none"
        }else{
            document.getElementsByClassName('ionex_file')[0].style.visibility="visible"
            document.getElementsByClassName('ionex_file')[0].style.display="flex"
        }
    })
     
}

const loader=()=>{
    // document.getElementById("submit").addEventListener('click',()=>{
    document.getElementById("submit").style.visibility="hidden"
    document.getElementById("submit").style.display="none"
    document.getElementsByClassName("lds-roller")[0].style.visibility="visible"
    document.getElementsByClassName("lds-roller")[0].style.display="flex"
    // })
}