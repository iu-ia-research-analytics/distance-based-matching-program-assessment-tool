##########Packages
library(shiny)
library(shinythemes)
library(DT)
library(tidyverse)
library(dplyr)
library(shinycssloaders)
library(MatchIt)
library(optmatch)
library(gridExtra)
library(clue)
library(plyr)
library(ggplot2)
library(packrat)
library(rsconnect)


##########Functions

#change data type function
change.data.type<-function(data,change.to.numeric,change.to.character){
  data<-data%>%
    mutate_at(change.to.numeric,as.numeric)%>%
    mutate_at(change.to.character,as.character)
  return(data)
}

#not in
`%not in%` <- Negate(`%in%`)


#normalize variable function (divide by max in absolute value for scaling)
normalize.max<-function(x){
  x/max(abs(x),na.rm=TRUE)
}

#function to create data frame of custom distances
custom.gower.df<-function(data,id,treatment,treatment.level,matching,matching.weight){
  data<-data%>%
    mutate(TREATMENT_INDICATOR=ifelse(data[,treatment]==treatment.level,'TREATMENT','CONTROL'))%>%
    mutate(POSITION=1:nrow(data))%>%
    mutate_at(matching[sapply(data[,matching],is.numeric)],normalize.max)
  data.treatment<-data%>%
    dplyr::filter(TREATMENT_INDICATOR=='TREATMENT')%>%
    select(id,matching)
  data.control<-data%>%
    dplyr::filter(TREATMENT_INDICATOR=='CONTROL')%>%
    select(id,matching)
  list.distances<-expand.grid(data.treatment[,id],
                              data.control[,id])
  names(list.distances)<-c('TREATMENT_ID','CONTROL_ID')
  list.distances$GOWER_SUM<-rep(0,nrow(list.distances))
  data_tmp<-data.frame()
  treatment_m<-data.frame()
  control_m<-data.frame()
  numeric_ind=NULL
  for(j in 2:ncol(data.treatment)){
    treatment_m<-data.treatment[,c(1,j)]
    control_m<-data.control[,c(1,j)]
    data_tmp<-list.distances[,c(1,2)]%>%
      inner_join(treatment_m,by=c('TREATMENT_ID'=id))%>%
      inner_join(control_m,by=c('CONTROL_ID'=id))
    data_tmp$VAR_TREATMENT<-data_tmp[,3]
    data_tmp$VAR_CONTROL<-data_tmp[,4]
    data_tmp<-data_tmp[,-c(3,4)]
    numeric_ind<-is.numeric(data_tmp$VAR_TREATMENT)
    if(numeric_ind==TRUE){
      data_tmp<-data_tmp%>%
        mutate(gower=if_else(!is.na(VAR_TREATMENT)&!is.na(VAR_CONTROL),
                             abs(VAR_TREATMENT-VAR_CONTROL),
                             1)*matching.weight[j-1])
    }
    if(numeric_ind==FALSE){
      data_tmp<-data_tmp%>%
        mutate(gower=if_else(!is.na(VAR_TREATMENT)&!is.na(VAR_CONTROL),
                             if_else(VAR_TREATMENT==VAR_CONTROL,
                                     0,
                                     1),
                             1)*matching.weight[j-1])
    }
    list.distances<-list.distances%>%
      inner_join(data_tmp,by=c("TREATMENT_ID"="TREATMENT_ID","CONTROL_ID"="CONTROL_ID"))%>%
      mutate(GOWER_SUM=GOWER_SUM+gower)%>%
      select(TREATMENT_ID,CONTROL_ID,GOWER_SUM)
  }
  list.distances<-list.distances%>%
    mutate(DISTANCE=GOWER_SUM/sum(matching.weight))%>%
    select(TREATMENT_ID,CONTROL_ID,DISTANCE)
  return(list.distances)
}

#function to create matrix of custom distances
custom.gower.matrix<-function(data,id,treatment,treatment.level,matching,matching.weight){
  data<-data%>%
    mutate(TREATMENT_INDICATOR=ifelse(data[,treatment]==treatment.level,'TREATMENT','CONTROL'))%>%
    mutate(POSITION=1:nrow(data))%>%
    mutate_at(matching[sapply(data[,matching],is.numeric)],normalize.max)
  data.treatment<-data%>%
    dplyr::filter(TREATMENT_INDICATOR=='TREATMENT')%>%
    select(id,matching)
  data.control<-data%>%
    dplyr::filter(TREATMENT_INDICATOR=='CONTROL')%>%
    select(id,matching)
  list.distances<-expand.grid(data.treatment[,id],
                              data.control[,id])
  names(list.distances)<-c('TREATMENT_ID','CONTROL_ID')
  list.distances$GOWER_SUM<-rep(0,nrow(list.distances))
  data_tmp<-data.frame()
  treatment_m<-data.frame()
  control_m<-data.frame()
  numeric_ind=NULL
  for(j in 2:ncol(data.treatment)){
    treatment_m<-data.treatment[,c(1,j)]
    control_m<-data.control[,c(1,j)]
    data_tmp<-list.distances[,c(1,2)]%>%
      inner_join(treatment_m,by=c('TREATMENT_ID'=id))%>%
      inner_join(control_m,by=c('CONTROL_ID'=id))
    data_tmp$VAR_TREATMENT<-data_tmp[,3]
    data_tmp$VAR_CONTROL<-data_tmp[,4]
    data_tmp<-data_tmp[,-c(3,4)]
    numeric_ind<-is.numeric(data_tmp$VAR_TREATMENT)
    if(numeric_ind==TRUE){
      data_tmp<-data_tmp%>%
        mutate(gower=if_else(!is.na(VAR_TREATMENT)&!is.na(VAR_CONTROL),
                             abs(VAR_TREATMENT-VAR_CONTROL)
                             ,
                             1)*matching.weight[j-1])
    }
    if(numeric_ind==FALSE){
      data_tmp<-data_tmp%>%
        mutate(gower=if_else(!is.na(VAR_TREATMENT)&!is.na(VAR_CONTROL),
                             if_else(VAR_TREATMENT==VAR_CONTROL,
                                     0,
                                     1),
                             1)*matching.weight[j-1])
    }
    list.distances<-list.distances%>%
      inner_join(data_tmp,by=c("TREATMENT_ID"="TREATMENT_ID","CONTROL_ID"="CONTROL_ID"))%>%
      mutate(GOWER_SUM=GOWER_SUM+gower)%>%
      select(TREATMENT_ID,CONTROL_ID,GOWER_SUM)
  }
  list.distances<-list.distances%>%
    mutate(DISTANCE=GOWER_SUM/sum(matching.weight))%>%
    select(TREATMENT_ID,CONTROL_ID,DISTANCE)
  id_position<-data%>%
    select(id,POSITION)
  list.distances.format<-list.distances%>%
    inner_join(id_position,by=c('TREATMENT_ID'=id))%>%
    dplyr::rename(TREATMENT_POSITION=POSITION)%>%
    inner_join(id_position,by=c("CONTROL_ID"=id))%>%
    dplyr::rename(CONTROL_POSITION=POSITION)%>%
    select(DISTANCE,TREATMENT_POSITION,CONTROL_POSITION)
  treatment_position<-data.frame(list.distances.format$TREATMENT_POSITION)%>%
    distinct()
  colnames(treatment_position)<-'TREATMENT_POSITION'
  control_position<-data.frame(list.distances.format$CONTROL_POSITION)%>%
    distinct()
  colnames(control_position)<-'CONTROL_POSITION'
  list.distances.format<-spread(list.distances.format,key=CONTROL_POSITION,value=DISTANCE)%>%
    select(-TREATMENT_POSITION)
  list.distances.matrix<-as.matrix(list.distances.format)
  rownames(list.distances.matrix)<-treatment_position[["TREATMENT_POSITION"]]
  return(list.distances.matrix)
}

#variable ratio matching function
varratio_matchit<-function(data,id,treatment,treatment.level,matching,weight,distance.matrix,method,avgratio,maxratio){
  data<-data%>%
    mutate(TREATMENT_INDICATOR=ifelse(data[,treatment]==treatment.level,1,0))
  m.out<-matchit(formula=TREATMENT_INDICATOR~1,
                 data=data, 
                 distance=distance.matrix,
                 method=method,
                 replace=FALSE, 
                 ratio=maxratio) 
  match_data<-match.data(m.out)
  distance_matrix_var_ratio<-custom.gower.matrix(
    data=match_data,
    id=id,
    treatment=treatment,
    treatment.level=treatment.level,
    matching=matching,
    matching.weight = weight
  )
  
  avg_matches<-avgratio
  max_matches<-maxratio
  min_matches<-1
  
  total.base.row<-nrow(distance_matrix_var_ratio)*max_matches
  total.col.sinks<-total.base.row-nrow(distance_matrix_var_ratio)*avg_matches
  total.col<-ncol(distance_matrix_var_ratio)+total.col.sinks
  total.row.sinks<-total.col-total.base.row
  
  matrix_tmp<-matrix(rep(NA,nrow(distance_matrix_var_ratio)*ncol(distance_matrix_var_ratio)),nrow=nrow(distance_matrix_var_ratio),byrow=TRUE)
  distance_matrix_assignment<-matrix(rep(NA,nrow(distance_matrix_var_ratio)*ncol(distance_matrix_var_ratio)),nrow=nrow(distance_matrix_var_ratio),byrow=TRUE)
  
  for(i in 1:(max_matches)){
    if(i==1){
      matrix_tmp<-distance_matrix_var_ratio
      rownames(matrix_tmp)<-paste(rownames(distance_matrix_var_ratio),'_',as.character(i),sep="")
      distance_matrix_assignment<-rbind(distance_matrix_assignment,matrix_tmp)
      distance_matrix_assignment<-distance_matrix_assignment[-c(seq(1,nrow(distance_matrix_var_ratio),by=1)),]
    }
    if(i!=1){
      matrix_tmp<-distance_matrix_var_ratio
      rownames(matrix_tmp)<-paste(rownames(distance_matrix_var_ratio),'_',as.character(i),sep="")
      distance_matrix_assignment<-rbind(distance_matrix_assignment,matrix_tmp)
    }
  }
  
  distance_matrix_assignment<-distance_matrix_assignment+50
  
  col.sinks.fudge1<-1000000+abs(rnorm(total.col.sinks*nrow(distance_matrix_var_ratio),50,10))
  col.sinks.fudge2<-0+abs(rnorm(total.col.sinks*nrow(distance_matrix_var_ratio)*(max_matches-1),5,5))
  
  col.sinks.values<-c(col.sinks.fudge1,col.sinks.fudge2)
  col.sinks.matrix<-matrix(col.sinks.values,nrow=total.base.row,ncol=total.col.sinks,byrow=TRUE)
  rownames(col.sinks.matrix)<-rownames(distance_matrix_assignment)
  colnames(col.sinks.matrix)<-paste("cs","_",c(1:ncol(col.sinks.matrix)),sep="")
  
  
  cost<-cbind(distance_matrix_assignment,col.sinks.matrix)
  
  rm(list=c('matrix_tmp',
            'distance_matrix_assignment',
            'col.sinks.values',
            'col.sinks.matrix'
  ))
  gc()
  
  assignment_solution<-solve_LSAP(cost,maximum=FALSE)
  
  assignment_vector<-
    data.frame(cbind(seq_along(assignment_solution), assignment_solution))
  assignment_vector$assignment_solution<-as.numeric(assignment_vector$assignment_solution)
  
  crosswalk<-data.frame(
    cbind(1:ncol(distance_matrix_var_ratio),
          colnames(distance_matrix_var_ratio))
  )
  names(crosswalk)<-c('COLUMN_SEQUENCE','ROWNUMBER')
  crosswalk$ROWNUMBER<-as.numeric(crosswalk$ROWNUMBER)
  crosswalk$COLUMN_SEQUENCE<-as.numeric(crosswalk$COLUMN_SEQUENCE)
  
  list1<-assignment_vector%>%
    dplyr::filter(assignment_solution<=total.base.row)%>%
    inner_join(crosswalk,by=c('assignment_solution'='COLUMN_SEQUENCE'))%>%
    select(ROWNUMBER)
  
  list2<-data.frame(as.numeric(rownames(distance_matrix_var_ratio)))
  names(list2)<-c('ROWNUMBER')
  
  list_use<-rbind(list1,list2)
  list_use$ROWNUMBER<-as.numeric(list_use$ROWNUMBER)
  
  data_balance<-
    match_data%>%
    mutate(ROWNUMBER=1:nrow(match_data))%>%
    inner_join(list_use,by=c('ROWNUMBER'='ROWNUMBER'))
  return(data_balance)
}

##########UI
ui<-navbarPage("Distance-Based Matching Program Assessment Tool",
               tabPanel("Load CSV File",
                        fileInput("file","Choose CSV File",accept=".csv"),
                        h5("How to Cite This Tool:"),
                        h6("Deom, Gina. Fiorini, Stefano. (2023). A Rapid Approach to Learning Analytics: A Distance-Based Program Assessment Tool. LAK 23 Accepted Paper. https://www.solaresearch.org/events/lak/lak23/accepted-papers/"),
               ),
               tabPanel("Variable Parameters",
                        fluidPage(theme=shinytheme("cerulean"),
                                  sidebarLayout(
                                    sidebarPanel(
                                      h2("Select Variables:"),
                                      selectInput(inputId="id",label="Unique Row Identifier:",choices=c('None Yet'),multiple=FALSE),
                                      selectInput(inputId="treatment",label="Treatment Variable:",choices=c('None Yet'),multiple=FALSE),
                                      selectInput(inputId="matching",label="Matching Variables:",choices=c('None Yet'),multiple=TRUE),
                                      selectInput(inputId="outcome",label="Outcome Variables:",choices=c('None Yet'),multiple=TRUE),
                                      selectInput(inputId="matching.change.numeric",label="Select Character Variables to Change to Numeric:",choices=c('None Yet'),multiple=TRUE),
                                      selectInput(inputId="matching.change.character",label="Select Numeric Variables to Change to Character:",choices=c('None Yet'),multiple=TRUE),
                                      submitButton("Submit", icon("arrows-rotate"))
                                    ),
                                    mainPanel(
                                      h4("Treatment Level & Data Types:"),
                                      uiOutput("treatment.levels"),
                                      submitButton("Submit", icon("arrows-rotate")),
                                      h4("Review Matching and Outcome Variable Data Types"),
                                      verbatimTextOutput("data.var.type")
                                    )
                                  )
                        )),
               tabPanel("Matching Parameters",
                        fluidPage(theme=shinytheme("cerulean")),
                        sidebarLayout(
                          sidebarPanel(
                            h2("Set Matching Variable Weights:"),
                            uiOutput("sliders"),
                            h2("Set Matching Parameters:"),
                            #selectInput(inputId="matching.type",label="Matching Type",choices=c('Fixed Ratio Matching'#,'Variable Ratio Matching (only use if available.ratio>=3)'
                            #                                                                    ),selected='Fixed Ratio Matching',multiple=FALSE),
                            selectInput(inputId="matching.algorithm",label="Matching Algorithm",choices=c('nearest','optimal'),selected='nearest',multiple=FALSE),
                            submitButton("Submit", icon("arrows-rotate"))
                          ),
                          mainPanel(
                            h4("Review Matching Variable Weights"),
                            verbatimTextOutput("weight.table"),
                            h4("Set Matching Ratios"),
                            verbatimTextOutput("treatment.control.nsize"),
                            uiOutput("ratiosliders"),
                            submitButton("Submit", icon("arrows-rotate")))
                        )),
               tabPanel("Matching Balance",
                        fluidPage(theme=shinytheme("cerulean")),
                        sidebarLayout(
                          sidebarPanel(
                            h2("Select Variable to Evaluate Balance"),
                            shinycssloaders::withSpinner(uiOutput("matching.var.balance")),
                            submitButton("Submit", icon("arrows-rotate")),
                            downloadButton("downloadparameters","Download Match Parameters"),
                            downloadButton("downloadmatchdata", "Download Matched Data"),
                            downloadButton("downloaddistances","Download Pairwise Distances*"),
                            h6("*Note: Please avoid clicking the download button multiple times. After clicking the button once, please allow for 1-2 minutes for the file to start appearing in downloads. Sometimes the file is large and can take some time to download."),
                            h6("*Note: After downloading files, move the files to a secure location or file server that is approved for the storage of applicable data at IU and delete the files from your download history.")
                          ),
                          mainPanel(
                            h2("Matching Summary"),
                            h4('Sample Sizes:'),
                            shinycssloaders::withSpinner(
                              verbatimTextOutput("match_sample_size")),
                            h4('Balance Before Matching (Press Submit):'),
                            shinycssloaders::withSpinner(plotOutput("plot_b4_match")),
                            h4('Balance After Matching (Press Submit):'),
                            shinycssloaders::withSpinner(plotOutput("plot_after_match")),
                            h4('Statistical Tests for Balance After Matching (Press Submit):'),
                            shinycssloaders::withSpinner(verbatimTextOutput("balance_test"))
                          )
                          
                        )
                        
               )
)

##########Server
server<-function(input,output,session){
  data_input <- reactive({
    file <- input$file
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(file$datapath, header = input$header)
  })
  observeEvent(input$file,{
    mytable<-read.csv(input$file$datapath)
    updateSelectInput(session,"id",label="Unique Row Identifier:",choices=c('Choose',colnames(mytable)))
    updateSelectInput(session,"treatment",label="Treatment Variable:",choices=c('Choose',colnames(mytable)))
    updateSelectInput(session,"matching",label="Matching Variables:",choices=colnames(mytable))
    updateSelectInput(session,"outcome",label="Outcome Variables:",choices=colnames(mytable))
    updateSelectInput(session,"matching.change.numeric",label="Select Character Variables to Change to Numeric:",choices=c('Choose',colnames(mytable[,unlist(lapply(mytable,is.character))])))
    updateSelectInput(session,"matching.change.character",label="Select Numeric Variables to Change to Character:",choices=c('Choose',colnames(mytable[,unlist(lapply(mytable,is.numeric))])))
  })
  output$treatment.levels<-renderUI({
    mytable<-read.csv(input$file$datapath) 
    selectInput(inputId="levels",label="Which level of the treatment variable identifies those who received the treatment?",choices=unique(mytable[[input$treatment]]),multiple=FALSE)
  })
  output$data.var.type<-renderPrint({
    mytable<-read.csv(input$file$datapath) 
    x.match<-change.data.type(data=mytable,
                              change.to.numeric=input$matching.change.numeric,
                              change.to.character=input$matching.change.character)
    x.match<-x.match%>%
      select(input$matching,input$outcome
      )
    sapply(x.match,typeof)
  })
  output$weights<-renderPrint({
    x.tmp<-change.data.type(data=data_input(),
                            change.to.numeric=input$matching.change.numeric,
                            change.to.character=input$matching.change.character)
    x.tmp<-x.tmp%>%
      select(input$matching)
    sapply(x.tmp,default.weights)
  })
  output$sliders<-renderUI({
    num.match.var<-length(input$matching)
    lapply(1:num.match.var,function(i){
      sliderInput(inputId = paste0("match", i), label = paste("Weight for ", input$matching[i]),
                  min = 1, max = 10, value = 1, step = 1)
    })
  })
  output$weight.vector<-renderPrint({
    num.match.var<-length(input$matching)
    sapply(1:num.match.var, function(i) {
      input[[paste0("match", i)]]
    })
  })
  output$weight.table<-renderPrint({
    num.match.var<-length(input$matching)
    match.var<-input$matching
    match.weights<-sapply(1:num.match.var, function(i) {
      input[[paste0("match", i)]]
    })
    data.frame(cbind(match.var,match.weights))
  })
  data_input_update<-reactive({
    mytable<-read.csv(input$file$datapath)
    x.tmp<-change.data.type(data=mytable,
                            change.to.numeric=input$matching.change.numeric,
                            change.to.character=input$matching.change.character)
    x.tmp$TREATMENT_INDICATOR=ifelse(x.tmp[,input$treatment]==input[["levels"]],1,0)
    x.tmp$TREATMENT_INDICATOR_LABEL=ifelse(x.tmp[,input$treatment]==input[["levels"]],'TREATMENT','CONTROL')
    x.tmp
  })
  data_input_treatment<-reactive({
    mytable<-read.csv(input$file$datapath)
    x.tmp<-change.data.type(data=mytable,
                            change.to.numeric=input$matching.change.numeric,
                            change.to.character=input$matching.change.character)
    x.tmp[x.tmp[[input$treatment]] %in% input[["levels"]], ]
  })
  data_input_control<-reactive({
    mytable<-read.csv(input$file$datapath)
    x.tmp<-change.data.type(data=mytable,
                            change.to.numeric=input$matching.change.numeric,
                            change.to.character=input$matching.change.character)
    x.tmp[x.tmp[[input$treatment]] %not in% input[["levels"]], ]
  })
  
  output$treatment.control.nsize<-renderPrint({
    n.treatment<-nrow(data_input_treatment())
    n.control<-nrow(data_input_control())
    available.ratio<-round(n.control/n.treatment,2)
    cbind(n.treatment,n.control,available.ratio)
  })
  output$ratiosliders<-renderUI({
    #req(input$matching.type)
    ntreatment<-nrow(data_input_treatment())
    ncontrol<-nrow(data_input_control())
    available_ratio<-round(ncontrol/ntreatment,2)
    ratio_max<-min(100,available_ratio)
    # if(input$matching.type=='Fixed Ratio Matching'){
    sliderInput(inputId="ratio.f",label='Select the number controls to be matched to each treatment.',min=1,max=ratio_max,value=1,step=1)
    # }
    # else if(input$matching.type=='Variable Ratio Matching (only use if available.ratio>=3)'){
    #   list(sliderInput(inputId="ratio.vra",label='Select the average number of controls to be matched to each treatment.',min=2,max=ratio_max,value=2,step=1),
    #     sliderInput(inputId="ratio.vra.max",label='Select the max number of controls that can be matched to any one treatment.',min=2,max=ratio_max,value=3,step=1))
    # }
  })
  output$matching.var.balance<-renderUI({
    selectInput(inputId="matching.filter",label="Select Matching Variable:",choices=unique(input$matching),selected=input$matching[1],multiple=FALSE)
  })
  matched_data<-reactive({
    # if(input$matching.type=='Fixed Ratio Matching'){
    data_tmp<-data_input_update()
    #data_tmp$TREATMENT_INDICATOR=ifelse(data_tmp[,input$treatment]==input[["levels"]],1,0)
    num.match.var<-length(input$matching)
    weight.vector<-sapply(1:num.match.var, function(i) {
      input[[paste0("match", i)]]
    })
    n_treatment<-nrow(data_tmp%>%dplyr::filter(TREATMENT_INDICATOR==1))
    distance_matrix<-custom.gower.matrix(
      data=data_tmp,
      id=input$id,
      treatment=input$treatment,
      treatment.level=input[["levels"]],
      matching=input$matching,
      matching.weight = weight.vector
    )
    m.out<-matchit(formula=TREATMENT_INDICATOR~1, #if using a distance matrix, do not need to supply list of variables here after ~ to match on
                   data=data_tmp, #data set
                   distance=distance_matrix, #distance matrix from step above
                   method=ifelse(n_treatment<=1,'optimal',input$matching.algorithm), #matching method (i.e. nearest or optimal)
                   replace=FALSE, #TRUE/FALSE identifying if matching done with replacement; i.e. if control units can be matched to more than one treatment (FALSE matching done without replacement)
                   ratio=input[["ratio.f"]]) #ratio of number of controls to number of treatments
    match.data(m.out)
    # }
    # else if(input$matching.type=='Variable Ratio Matching (only use if available.ratio>=3)'){
    #   data_tmp<-data_input_update()
    #   num.match.var<-length(input$matching)
    #   weight.vector<-sapply(1:num.match.var, function(i) {
    #     input[[paste0("match", i)]]
    #   })
    #   distance_matrix<-custom.gower.matrix(
    #     data=data_tmp,
    #     id=input$id,
    #     treatment=input$treatment,
    #     treatment.level=input[["levels"]],
    #     matching=input$matching,
    #     matching.weight = weight.vector
    #   )
    #   varratio_matchit(data=data_tmp,
    #                    id=input$id,
    #                    treatment=input$treatment,
    #                    treatment.level=input[["levels"]],
    #                    matching=input$matching,
    #                    weight=weight.vector,
    #                    distance.matrix=distance_matrix,
    #                    method=input$matching.algorithm,
    #                    avgratio=input[["ratio.vra"]],
    #                    maxratio=input[["ratio.vra.max"]]
    #                    
    #   )
    #  }
  })
  distance_vector<-reactive({
    data_tmp<-data_input_update()
    num.match.var<-length(input$matching)
    weight.vector<-sapply(1:num.match.var, function(i) {
      input[[paste0("match", i)]]
    })
    custom.gower.df(
      data=data_tmp,
      id=input$id,
      treatment=input$treatment,
      treatment.level=input[["levels"]],
      matching=input$matching,
      matching.weight = weight.vector
    )
  })
  match_filter_selection<-reactive({
    input$matching.filter
  })
  matched_data_print<-renderPrint({
    head(matched_data())
  })
  data_tmp<-reactive({
    tmp<-data_input_update()
    tmp_variable<-ifelse(is.na(match_filter_selection()),input$matching[1],match_filter_selection())
    tmp<-data.frame(tmp[,c(tmp_variable,'TREATMENT_INDICATOR_LABEL')])
    names(tmp)<-c('VARIABLE','TREATMENT')
    tmp
  })
  data_tmp_after<-reactive({
    tmp<-matched_data()
    tmp_variable<-ifelse(is.na(match_filter_selection()),input$matching[1],match_filter_selection())
    tmp<-data.frame(tmp[,c(tmp_variable,'TREATMENT_INDICATOR_LABEL')])
    names(tmp)<-c('VARIABLE','TREATMENT')
    tmp
  })
  numeric_flag<-reactive({
    data_tmp<-data_input_update()
    tmp_variable<-ifelse(is.na(match_filter_selection()),input$matching[1],match_filter_selection())
    tmp<-data_tmp%>%
      select(tmp_variable)
    #select(input$matching.filter)
    is.numeric(tmp[1,1])
  })
  output$match_sample_size<-renderPrint({
    data_tmp<-matched_data()
    #data_tmp<-data_input_update()
    n.treatment<-nrow(data_tmp[data_tmp[[input$treatment]] %in% input[["levels"]], ])
    #data_tmp
    n.control<-nrow(data_tmp[data_tmp[[input$treatment]] %not in% input[["levels"]], ])
    #data_tmp
    ratio<-n.control/n.treatment
    cbind(n.treatment,n.control,ratio)
  })
  output$plot_b4_match<-renderPlot({
    if(numeric_flag()==FALSE){
      ggplot(data_tmp(),aes(x=VARIABLE,group=TREATMENT))+
        geom_bar(aes(y=..prop..,fill=factor(..x..)
        ),stat="count")+
        geom_text(aes(label=scales::percent(..prop..),
                      y=..prop..),stat="count",vjust=-.5)+
        labs(y="Percent")+
        facet_grid(~TREATMENT)+
        scale_y_continuous(labels=scales::percent)+
        theme(legend.position="none")
    }
    else if(numeric_flag()==TRUE){
      ggplot(data_tmp(),aes(x=TREATMENT,y=VARIABLE,fill=TREATMENT))+
        geom_boxplot()+
        theme(legend.position="none")
    }
  })
  output$plot_after_match<-renderPlot({
    if(numeric_flag()==FALSE){
      ggplot(data_tmp_after(),aes(x=VARIABLE,group=TREATMENT))+
        geom_bar(aes(y=..prop..,fill=factor(..x..)
        ),stat="count")+
        geom_text(aes(label=scales::percent(..prop..),
                      y=..prop..),stat="count",vjust=-.5)+
        labs(y="Percent")+
        facet_grid(~TREATMENT)+
        scale_y_continuous(labels=scales::percent)+
        theme(legend.position="none")
    }
    else if(numeric_flag()==TRUE){
      ggplot(data_tmp_after(),aes(x=TREATMENT,y=VARIABLE,fill=TREATMENT))+
        geom_boxplot()+
        theme(legend.position="none")
    }
  })
  output$balance_test<-renderPrint({
    if(!is.null(match_filter_selection())==TRUE){
      tmp<-data_tmp_after()
      if(numeric_flag()==FALSE){
        try(chisq.test(tmp$TREATMENT,tmp$VARIABLE),silent=TRUE)
      }
      else if(numeric_flag()==TRUE){
        try(t.test(tmp$VARIABLE~tmp$TREATMENT),silent=TRUE)
      }
    }
    else ""
  })
  matching_parameters<-reactive({
    num.match.var<-length(input$matching)
    num.outcome.var<-length(input$outcome)
    weight.vector<-sapply(1:num.match.var, function(i) {
      input[[paste0("match", i)]]
    })
    parameter<-c('Unique Row Identifier',
                 'Treatment Variable',
                 rep('Matching Variable',num.match.var),
                 rep('Outcome Variable',num.outcome.var),
                 'Matching Algorithm',
                 'Matching Ratio')
    variable<-c(input$id,
                input$treatment,
                input$matching,
                input$outcome,
                input$matching.algorithm,
                input[["ratio.f"]])
    weight_level<-c(NA,
                    input[["levels"]],
                    weight.vector,
                    NA,
                    NA,
                    NA)
    data_tmp<-matched_data()%>%
      select(c(input$matching,input$outcome))
    data_type_tmp<-data.frame(sapply(data_tmp,typeof))
    names(data_type_tmp)<-c('type')
    data_type_tmp<-data_type_tmp$type
    data_type<-c(NA,
                 NA,
                 data_type_tmp,
                 NA,
                 NA)
    parameter_summary<-data.frame(cbind(parameter,
                                        variable,
                                        weight_level,
                                        data_type))
    names(parameter_summary)<-c('PARAMETER',
                                'VALUE',
                                'WEIGHT_OR_LEVEL',
                                'DATA_TYPE')
    parameter_summary
  })
  output$matching_parameters_print<-renderPrint({
    head(matching_parameters())
  })
  output$downloadmatchdata<-downloadHandler(
    filename = function() {
      paste('Matched_Data', Sys.Date(), '.csv',sep='')
    },
    content = function(file) {
      write.csv(matched_data(),file,row.names=FALSE)
    },
    contentType="text/csv"
  )
  output$downloaddistances<-downloadHandler(
    filename = function() {
      paste('Distances', Sys.Date(), '.csv',sep='')
    },
    content = function(file) {
      write.csv(distance_vector(), file,row.names=FALSE)
    },
    contentType="text/csv"
  )
  output$downloadparameters<-downloadHandler(
    filename = function() {
      paste('Parameters', Sys.Date(), '.csv',sep='')
    },
    content = function(file) {
      write.csv(matching_parameters(), file,row.names=FALSE)
    },
    contentType="text/csv" 
  )
}

##########Create App
shinyApp(ui=ui,server=server)