#' Plot maps of indexes
#'
#' @param abundance_grid spatial table of abundance
#' @param biomass_grid spatial table of biomass
#' @param meanWEIGHT_grid spatial table of mean weight
#' @param map_range range of coordinates for the map
#' @param stratum reference stratum range (allowed values: "10,200","10,800","200,800","5,35","5,45")
#' @importFrom terra wrap unwrap
#' @export
plotGFCM <- function (abundance_grid,biomass_grid,meanWEIGHT_grid, map_range, stratum){

    if(FALSE) {
        map_range <- c(15,21,39.5,42.5)
    }

    sspp <- unique(abundance_grid$sspp)[1]
    names(map_range) <- c("xmin","xmax","ymin","ymax")
    summary.grid <- abundance_grid
    # metaDB <- mTATB

            x1 <- min(map_range[1])
            x2 <- max(map_range[2])
            y1 <- min(map_range[3])
            y2 <- max(map_range[4])

    #----------------------------------------
    # Saving maps of the ABUNDANCE index
    #----------------------------------------                                                                          # [summary.grid$meanNkm2>5]
    cat_lim <- as.numeric(quantile(summary.grid[!is.na(summary.grid$meanNkm2), "meanNkm2"], probs = seq(0, 1,0.20) ) )
    cate_label <-   round(cat_lim,1)
    cate_label_min <- cate_label
    cate_label_max <- cate_label[2:length(cate_label)]


    cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
    cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")
    cate_label <- cate_label[1: (length(cate_label)-1)]

    palette_cate <- heat.colors(length(cate_label))
    palette_cate <- palette_cate[length(palette_cate):1]

    cgpmgrid <- unwrap(cgpmgrid)
    continent <- unwrap(continent)
    cgpmgrid$abundance <- NA
    cgpmgrid$color <- "#FFFFFF"
    i=2
    for (i in 1:length(summary.grid[,1])){
        cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"meanNkm2"]
        if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
        }

        if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3] ) & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
        }

        if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
        }

        if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
        }

        if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
        }
    }


    lx =map_range["xmax"] - map_range["xmin"]
    ly =map_range["ymax"] - map_range["ymin"]
    ratio <- ly/lx
    img_width <- 30
    img_height <- img_width*ratio

    par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
    plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]) ,  xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste(paste("Abundance  (n/km^2)", sep=" "), sspp), asp=1)
    plot(cgpmgrid,col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")                                                                                 #data = centroidi_coords,
    if (stratum == "200,800" | stratum == 3){
        plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
    }
    plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)
    legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Abundance (n/km^2)\n(classes based on percentiles)", bty = "n")
    # legend(map_range[1,"xmin"],map_range[1,"ymax"], cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Abundance (n/km^2)\n(classes based on percentiles)", bty = "n")
    box()
    p1 <- recordPlot()

    jpeg(file=paste(wd, "/output/",sspp," - GFCM GRID ABUNDANCE.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
     replayPlot(p1)
    dev.off()



    #----------------------------------------
    # Saving maps of the inverse of CV of ABUNDANCE index
    #----------------------------------------
    cat_lim <- quantile(summary.grid$inverse_cvNkm2[!is.na(summary.grid$inverse_cvNkm2)], probs = seq(0, 1,0.20) )
    cate_label <-  round(cat_lim,3)

    cate_label_min <- cate_label
    cate_label_max <- cate_label[2:length(cate_label)]


    cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
    cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")

    #cate_label[length(cate_label)] <- paste(">", cate_label_max[length(cate_label_max)] )
    cate_label <- cate_label[1: (length(cate_label)-1)]

    palette_cate <- heat.colors(length(cate_label))
    palette_cate <- palette_cate[length(palette_cate):1]

    cgpmgrid$invCV <- NA
    cgpmgrid$color <- "#FFFFFF"

    i=3
    for (i in 1:length(summary.grid[,1])){
        cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"inverse_cvNkm2"]

        if(all(is.na(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]])) ){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- "#FFFFFF"
        } else {

            if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
            }

            if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
            }

            if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
            }

            if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
            }

            if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
            }
        }
    }

    lx =map_range["xmax"] - map_range["xmin"]
    ly =map_range["ymax"] - map_range["ymin"]
    ratio <- ly/lx
    img_width <- 30
    img_height <- img_width*ratio

    par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
    plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste("Inverse of Abundance CV -", sspp) , asp=1)
    plot(cgpmgrid,col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")                                                                                 #data = centroidi_coords,
    if (stratum == "200,800" | stratum == 3){
        plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
    }
    plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)
    legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Inverse of Abundance CV\n(classes based on percentiles)", bty = "n")
    box()
    p2 <- recordPlot()

    jpeg(file=paste(wd, "/output/",sspp," - GFCM GRID ABUNDANCE Inverse CV.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
     replayPlot(p2)
    dev.off()



    ###################
    ###################
    ###             ###
    ###   BIOMASS   ###
    ###             ###
    ###################
    ###################


    summary.grid <- biomass_grid
    # metaDB <- mTATB

    #----------------------------------------
    # Saving maps of the BIOMASS index
    #----------------------------------------

    cat_lim <- quantile(summary.grid[!is.na(summary.grid$meanBkm2), "meanBkm2"], probs = seq(0, 1,0.20) )
    cate_label <-  round(cat_lim,3)

    cate_label_min <- cate_label
    # cate_label_min[1] <- 0
    cate_label_max <- cate_label[2:length(cate_label)]


    cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
    cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")

    cate_label <- cate_label[1: (length(cate_label)-1)]
    # cate_label[length(cate_label)] <- paste(">", cate_label_max[length(cate_label_max)] )


    palette_cate <- heat.colors(length(cate_label))
    palette_cate <- palette_cate[length(palette_cate):1]


    cgpmgrid$meanBkm2 <- NA
    cgpmgrid$color <- "#FFFFFF"

    i=3
    for (i in 1:length(summary.grid[,1])){
        cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"meanBkm2"]

        if(all(is.na(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]])) ){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- "#FFFFFF"
        } else {

            if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2] ) & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
            }

            if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
            }

            if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
            }

            if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
            }

            if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
            }
        }
    }


    lx =map_range["xmax"] - map_range["xmin"]
    ly =map_range["ymax"] - map_range["ymin"]
    ratio <- ly/lx
    img_width <- 30
    img_height <- img_width*ratio

    par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
    plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste("Biomass (kg/km^2) -", sspp), asp=1 )
    plot(cgpmgrid,col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")                                                                                 #data = centroidi_coords,
    if (stratum == "200,800" | stratum == 3){
        plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
    }
    plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)

    legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Biomass (kg/km^2)\n(classes based on percentiles)", bty = "n")
    box()
    p3 <- recordPlot()

    jpeg(file=paste(wd, "/output/",sspp," - GFCM GRID BIOMASS.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
replayPlot(p3)
    dev.off()




    ###################
    ###################
    ###             ###
    ###     MIW     ###
    ###             ###
    ###################
    ###################

    summary.grid <- meanWEIGHT_grid
    # metaDB <- mTATB

    #----------------------------------------
    # Saving maps of the MIW
    #----------------------------------------                                                                          # [summary.grid$meanNkm2>5]
    cat_lim <- as.numeric(quantile(summary.grid[!is.na(summary.grid$MIW), "MIW"], probs = seq(0, 1,0.20) ) )
    #cat_lim[1] <- 0
    #  cate_label[len] <- ifelse( round(cat_lim[len],0) != 0, round(cat_lim[len],0), ifelse( round(cat_lim[len],1) != 0,  round(cat_lim[len],1),  round(cat_lim[len],2)))
    cate_label <-   round(cat_lim,3)
    cate_label_min <- cate_label
    #cate_label_min[1] <- 0
    cate_label_max <- cate_label[2:length(cate_label)]


    cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
    cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")

    #cate_label[length(cate_label)] <- paste(">", cate_label_max[length(cate_label_max)] )
    cate_label <- cate_label[1: (length(cate_label)-1)]

    palette_cate <- heat.colors(length(cate_label))
    palette_cate <- palette_cate[length(palette_cate):1]

    cgpmgrid$MIW <- NA
    cgpmgrid$color <- "#FFFFFF"

    i=3
    for (i in 1:length(summary.grid[,1])){
        cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"MIW"]

        if(all(is.na(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]])) ){
            cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- "#FFFFFF"
        } else {

            if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
            }

            if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
            }

            if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
            }

            if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
            }

            if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
                cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
            }
        }
    }


    lx =map_range["xmax"] - map_range["xmin"]
    ly =map_range["ymax"] - map_range["ymin"]
    ratio <- ly/lx
    img_width <- 30
    img_height <- img_width*ratio


    par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
    plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]) ,  xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste(paste("Mean Individual Weight  (kg)", sep=" "), sspp), asp=1)
    plot(cgpmgrid, col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
    if (stratum == "200,800" | stratum == 3){
        plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
    }
    plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)

    legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "MIW (kg)\n(classes based on percentiles)", bty = "n")

    box()
    p4 <- recordPlot()



    jpeg(file=paste(wd, "/output/",sspp," - GFCM GRID MIW.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
replayPlot(p4)
    dev.off()


}
