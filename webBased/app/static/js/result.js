const handleFig=(prob)=>{
    console.log(JSON.parse(prob))
}

handleFig('{{ fig | tojson | safe}}');