#FinalProject
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(htmltools)
library(leaflet)
library(rgdal)
library(geojsonio)
library(jsonlite)
library(dplyr)
#source("jaxmat.R")   #for displaying mathematics
stylesheet <- tags$head(tags$style(HTML('
    .main-header .logo {
      font-family: "Georgia", Times, "Times New Roman", serif;
      font-weight: bold;
      font-size: 24px;
    }
  ')
))
  
#The user interface
header <- dashboardHeader(title = "COVID Deaths by Region",
                          titleWidth = 500)
sidebar <- dashboardSidebar(disable = TRUE)


body <- dashboardBody(
  fluidRow(stylesheet,
    column(width=1,
      actionBttn("btnq","Click me!"),
      h1(""),
      actionBttn("btngraph","GRAPH"),
      h3(""),
      actionBttn("btnmap","MAP"),
      h3(""),
      actionBttn("btncircles","CIRCLES"),
      h3(""),
      actionBttn("btntable","TABLE"),
      h3(""),
      selectInput("regions", "Regions", NULL),
    ),
    column(width = 6,
      leafletOutput("mymap"),
      h3("Cumulative COVID deaths After 200 Days"),
      plotOutput("cases"),
    ),
    column(width = 4,
           leafletOutput("myusmap"),
           h3("Confirmed Cases and Deaths by State"),
           numericInput("maxrows", "Rows to show", 15),
           verbatimTextOutput("rawtable"),
           
           
           
    ),
  )
)
ui <- dashboardPage(header, sidebar, body, skin = "blue") #other colors available

COVID <- read.csv("owid-covid-data.csv")
nRegion <- nrow(COVID)/150; nRegion
ndays <- 200
COVID[1:ndays,4]
#Functions to get access to the numeric data for a region
extract <- function(x) (COVID[(1+84*(x-1)):(ndays+84*(x-1)),4])
avg <- function(x) mean(extract(x))            #their average
name <- function(x) paste(COVID[1+84*(x-1),1]," ",COVID[1+84*(x-1),2])
extract(33)
name(33)
A <- matrix(ndays*nRegion,nrow = ndays, ncol = nRegion)      #make an empty matrix
#Center and rescale the data, one column per region
for (i in 1:nRegion){
  A[,i] <- (extract(i)-avg(i))/sqrt(ndays)
}
#Make a symmetric square matrix
S <- A%*%t(A);
Eig <- eigen(S)
Eig$values  #first few Coefficients are much larger than the others
P <- -Eig$vectors  #change the sign to get more positive coefficients



#Now we can express each column of A relative to the new basis of eigenvectors. 
#The complete set of coefficients will provide perfect reconstruction
PInv <- solve(P)
A.eig <- PInv%*%A
A.eig[,33]

#We can reconstruct the original data for one region using all the eigenvectors
k <- 33
name(k)
recon <- as.numeric(P%*%A.eig*sqrt(ndays)+avg(k))
#plot(1:ndays,recon,xlab = name(k), type = "l")
#points(1:ndays,extract(k))


server <- function(session, input, output) {
  cityDF <- sph.makeCityDF("covid.csv")
  usaDF <- sph.makeCityDF("usa.csv")
  updateSelectInput(session, "selstart", choices = cityDF$Name)
  updateSelectInput(session, "seldest", choices = cityDF$Name)
  colors <- c("black", "red", "green", "blue", "magenta", "gray", "orange", "brown", "purple", "cyan")
  color <- colorFactor(c("yellow", "orange", "red"), domain = c("7", "5", "3"))
  colorUS <- colorFactor(c("yellow", "orange", "purple", "red"), domain = c("7", "4", "3", "2"))
  index <- 1
  P <- -Eig$vectors  #change the sign to get more positive coefficients
  PInv <- solve(P)
  A.eig <- PInv%*%A
  regions <- character(nRegion)
  for (i in 1:nRegion)
    regions[i] <- trimws(name(i))
  updateSelectInput(session,"regions",choices = regions, selected = "Italy")
  D <- round(Eig$values)
  output$eigen <- renderPlot({
    plot(1:ndays, P[,1], ylim = c(-0.4,0.4),type = "l")
    for (i in 2:input$numshow)
      points(1:ndays, P[,i], type = "l",col = colors[i])
    
  })

 
  
  showCases <- function(casedata, rName){
    plot(1:ndays,casedata,xlab = rName, type = "l")
  }
  observeEvent(input$regions, {
    index <<- which.max(regions == input$regions)
    output$cases <- renderPlot(showCases(extract(index),input$regions))
    
  })
  
  output$rawtable <- renderPrint({
    orig <- options(width = 1000)
    print(tail(usaDF %>% select(c(Name, Confirmed, Deaths)), input$maxrows), row.names = FALSE)
    options(orig)
  })
  
  countries <- readOGR("time.geojson")
  qpal <- colorBin("Oranges", countries$current_deaths)
  coronavirus <- merge(cityDF, countries, by.x = "Name", by.y = "name")
  labels <- c("0-3", "3-7", "7+")
  drawMap <- function() {
    leaflet(data = cityDF) %>% 
      addTiles() %>%
      setView(lng = -43, lat = 0, zoom = 2) %>%
      addPolygons(data=countries, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 1,
                  color = ~qpal(current_deaths)) %>%
      addLegend(data=countries, "bottomleft", pal = qpal, values = ~current_deaths, title = "Total COVID Deaths",opacity = 1) %>%
      addCircleMarkers(lng=~Long,radius = ~as.numeric(sqrt(Deaths))/10, stroke=FALSE, # Circle stroke
                       fillOpacity=0.5, labelOptions = labelOptions(noHide = F, textsize = "15px"),color=~color(Mortality),
                       lat=~Lat, popup = paste("Confirmed Cases:", cityDF$Confirmed,sep="<br/>", "Deaths:", cityDF$Deaths),  label =  cityDF$Name)%>%
      addLegend("bottomright", pal = color, values = ~(Mortality), opacity = 2, title="Mortality Rate", labFormat = function(type, cuts, p) 
      {paste0(labels)}) 
    
    
    
  }
  output$mymap <- renderLeaflet ({
    drawMap()
  })
  
  
  labelUS <- c("0-2", "2-4", "4-6", "6-8")
  drawUSMap <- function() {
    leaflet(data = usaDF) %>% 
      addTiles() %>%
      setView(lng = -95, lat = 30, zoom = 3) %>%
      addCircleMarkers( lng=~Long,radius = ~sqrt(Deaths)/10, stroke=FALSE, # Circle stroke
                        fillOpacity=0.7, labelOptions = labelOptions(noHide = F, textsize = "15px"),color=~colorUS(Mortality),
                        lat=~Lat, popup = paste("Confirmed Cases:", usaDF$Confirmed,sep="<br/>", "Deaths:", usaDF$Deaths),  label =  usaDF$Name)%>%
      addLegend("bottomright", pal = colorUS, values = ~(Mortality), opacity = 2, title="Mortality Rate", labFormat = function(type, cuts, p) {
        paste0(labelUS)
      })
  }
  output$myusmap <- renderLeaflet ({
    drawUSMap()
  })
  
  observeEvent(input$btnq, {
    showModal(modalDialog(id = "yesno",
                          h3("Welcome to my app! Please click on the buttons below to learn how to interact with the map"),
                          footer = tagList(
                            modalButton("Close"),   #dismisses the dialog
                          )
    ))
  })
  observeEvent(input$btngraph, {
    showModal(modalDialog(id = "yesno",
                          h3(" Select the region from the dropdown to view a graph of cumulative COVID deaths over a 200 day period. "),
                          footer = tagList(
                            modalButton("Close"),   #dismisses the dialog
                          )
    ))
  })
  observeEvent(input$btncircles, {
    showModal(modalDialog(id = "yesno",
                          h3("1- Hover over any of the circles on the map to view the country name"),
                          h3("2- Clicking on the circle will show you more details on the COVID cases of the country"),
                          h3("Please note: the circles are different colors based on the Mortality Rate of the country (as indicated in the bottom right legend). They are also different sizes based on the number of COVID cases."),
                          footer = tagList(
                            modalButton("Close"),   #dismisses the dialog
                          )
    ))
  })
  
  observeEvent(input$btnmap, {
    showModal(modalDialog(id = "yesno",
                          h3("The countries on the world map are a different color based on the COVID deaths in the country (as indicated by the legend on the bottom left)"),
                          footer = tagList(
                            modalButton("Close"),   #dismisses the dialog
                          )
  
    ))
  })
  
  observeEvent(input$btntable, {
    showModal(modalDialog(id = "yesno",
                          h3("To change the number of states listed in the table, please edit the number under Rows to show "),
                          footer = tagList(
                            modalButton("Close"),   #dismisses the dialog
                          )
                          
    ))
  })
  
}

#Run the app
shinyApp(ui = ui, server = server)