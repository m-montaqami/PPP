const Btns=["btn1","btn2","btn3"]
const Items=["desc1","desc2","desc3"]
const descCtl=(prob1,prob2)=>{
        Items.map((e)=>{
            if(e==prob2){
                // console.log(e)
                
                document.getElementById(e).style.display="block";
                document.getElementById(e).style.visibility="visible";

            }else{
                document.getElementById(e).style.display="none";
                document.getElementById(e).style.visibility="hidden";
            }
        })    

        Btns.map((e)=>{
            if(e==prob1){
                document.getElementById(e).style.background="rgb(63, 82, 63)";

            }else{
                document.getElementById(e).style.background="rgb(146, 161, 146)";
            }
        })    

}