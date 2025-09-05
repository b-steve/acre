
#' Plot the detection function
#'
#' @param fit an object generated from the model fitting function [fit.acre] or
#'            the bootstrap process [boot.acre].
#' @param newdata data.frame; contains any covariates that will be used for all extended parameters (if not be skipped)
#' @param skip.extend.param character; skip extended parameter, for skipped extended parameters,
#'                          use its intercept as the value for this parameter
#' @param xlim a numeric vector with two elements as the range of x-axis.
#' @param ylim a numeric vector with two elements as the range of y-axis.
#' @param main a string as the main title of the plot.
#' @param xlab a string as the sub-title of x-axis.
#' @param ylab a string as the sub-title of y-axis.
#' @param col a string or a numeric vector indicates the color of the plotted line.
#' @param add a logical value indicates whether to add the lines into the existing plot.
#' @param ... 
#'
#' @return
#' @export
show.detfn <- function(fit, newdata = NULL, skip.extend.param = NULL, xlim = NULL, ylim = NULL,
                       main = NULL, xlab = NULL, ylab = NULL, col = NULL, add = FALSE, ...){

  det_fn = get_detfn(fit)
  
  if(is.null(xlab)) xlab = "Distance"
  if(is.null(ylab)) {
    if(det_fn == 'ss') ylab = "E(ss)"
    if(det_fn != 'ss') ylab = "Detection probability"
  }
  
  #get some essential component from output
  param_name = get_param_og(fit)
  param_name = param_name[-which(param_name == 'D')]
  data_param = get_data_param(fit)
  par_extend_name = get_par_extend_name(fit)
  param_values = get_coef(fit)
  ss_link = get_ss_link(fit)
  
  #deal with xlim and the values will be used as x-axis (distances)
  if (is.null(xlim)){
    buffer <- get_buffer(fit)
    #get_buffer return a vector with length of n.sessions
    #and if the buffer for on session is not provided, corresponding element is zero
    if (!any(buffer == 0)){
      x.max <- max(buffer)
    } else {
      x.max <- max(get_dist_theta(fit)$dx)
    }
    xlim <- c(0, x.max)
  }
  dists <- seq(xlim[1], xlim[2], length.out = 1000)
  
  #deal with skipped extended parameters and newdata
  
  if(is.null(newdata) & !is.null(skip.extend.param)){
    if(!all(skip.extend.param == par_extend_name[-which(par_extend_name == 'D')])){
      warning('No new data of covariates provided, all extended parameters will be skipped.')
    }
  }
  
  #if no newdata assigned, skip all extended parameters
  if(is.null(newdata) & is.null(skip.extend.param)){
    if(length(par_extend_name) > 0){
      message('No new data of covariates provided, only the intercept will be used for all extended parameters.')
      skip.extend.param = par_extend_name
    }
  }
  
  
  #set the values of every each detection function parameters
  det_param_input = vector('list', length(param_name))
  names(det_param_input) = param_name
  
  for(i in param_name){
    #get the basic information for this parameter
    tem = data_param[which(data_param$par == i),]
    link = tem$link
    n_col_full = tem$n_col_full
    n_col_mask = tem$n_col_mask
    
    values = param_values[[i]]
    
    if((i %in% par_extend_name) & (!i %in% skip.extend.param)){
      gam = get_gam(fit, i)
      det_param_input[[i]] = get_extended_par_value(gam, n_col_full, n_col_mask, values, newdata)
    } else {
      det_param_input[[i]] = values[1]
    }
    
    #use 'link' to back-transform the parameter's value
    
    det_param_input[[i]] = unlink.fun(link = link, value = det_param_input[[i]])
    
    names(det_param_input[[i]]) = NULL
  }
  
  #be careful, each component in det_param_input may have different length
  #but as.data.frame() will automatically solve this
  tem = as.data.frame(det_param_input)
  n_lines = nrow(tem)
  for(i in param_name) det_param_input[[i]] = tem[, i]
  
  probs = vector('list', n_lines)
  tem_det_par = vector('list', length(param_name))
  names(tem_det_par) = param_name
  
  for(i in 1:n_lines){
    #extract i'th value for each parameter from det_param_input
    for(j in param_name) tem_det_par[[j]] = det_param_input[[j]][i]
    probs[[i]] = det_prob(det_fn, tem_det_par, dists, ss_link)
  }
  
  #based on all probs values, determine ylim
  if(is.null(ylim)){
    if(det_fn != 'ss'){
      ylim = c(0, 1)
    } else {
      tem = do.call('c', probs)
      ylim = range(tem)
    }
  }
  
  if(is.null(col)){
    col = 1:n_lines
  } else {
    col = rep(col, length = n_lines)
  }
  
  if (!add){
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = 'i', yaxs = "r")
    title(main = "Detection function")
    axis(1)
    axis(2)
    box()
    if(det_fn != 'ss') abline(h = c(0, 1), col = "lightgrey")
    title(main = main, xlab = xlab, ylab = ylab)
  }
  
  for(i in 1:n_lines){
    lines(dists, probs[[i]], col = col[i], ...)
  }
  
}





#' Plotting an estimated density surface
#' 
#' Plots density surface estimated by a model fitted with the function [fit.acre]
#'
#' @param fit an object generated from the model fitting function "fit.acre()" or
#'            the bootstrap process "boot.acre()".
#' @param session The session with the detector array and invidual(s)
#'     to be plotted. Ignored if the `newdata` argument is
#'     provided.
#' @param new.data A data frame including new mask points and covariate
#'     values, from which to estimate and plot density estimates
#'     for. This allows, for example, estimates to be provided for new
#'     regions not included in the mask used to fit the model. Two
#'     columns, named `x` and `y`, must be included,
#'     providing the x- and y-coordinates of the new mask
#'     points. Additional columns must provide the covariates used to
#'     fit the model.
#' @param show.cv Logical. If true, the CV of the density estimate is
#'        plotted rather than the estimate itself. At present, this will
#'        only work if `newdata` is also provided.
#' @param D.cov 
#' @param x.pixels 
#' @param y.pixels 
#' @param zlim A numeric vector of length 2, giving the range of the density contours
#' @param scale By default, density is in animals per hectare. The 
#'        plotted values are multiplied by this argument, allowing 
#'        for user-specified units. For example, setting `scale = 100` 
#'        results in densities plotted as animals per square kilometre.
#' @param plot.contours Logical, if `TRUE`, contours are plotted. 
#' @param add a logical value indicates whether to add the lines into the existing plot
#' @param arg.col A numeric value, indicating the number of levels to stretch the color over
#' @param trap.plot 
#' @param ... 
#' @inheritParams show.detfn
#' @inheritParams read.acre
#'
#' @return
#' @export
#'
#' @examples
show.Dsurf <- function(fit, session = NULL, show.cv = FALSE, new.data = NULL, D.cov = NULL, xlim = NULL, ylim = NULL,
                        x.pixels = 50, y.pixels = 50, zlim = NULL, scale = 1, plot.contours = FALSE,
                        add = FALSE, convert.loc2mask= NULL, arg.col = 100, trap.plot = NULL, ...){
  
  pred = predict_D_for_plot(fit, session_select = ifelse(is.null(session), 1, session), 
                            new_data = new.data, D_cov = D.cov, xlim = xlim, ylim = ylim,
                            x_pixels = x.pixels, y_pixels = y.pixels, se_fit = show.cv,
                            convert.loc2mask= convert.loc2mask)
  #browser()
  mask = as.matrix(pred[, c('x', 'y')])
  if(!show.cv){
    D.mask = pred[,'est']
  } else {
    D.mask = pred[,'std'] / pred[,'est']
  }
  
  if (is.null(xlim)){
    xlim <- range(mask[, 'x'])
  }
  if (is.null(ylim)){
    ylim <- range(mask[, 'y'])
  }
  

  unique.x <- sort(unique(mask[, 1]))
  unique.y <- sort(unique(mask[, 2]))
  z <- squarify(mask, D.mask)
  

  if(!show.cv){
    z <- scale*z
  }
  
  if (is.null(zlim)){
    zlim <- c(0, max(z, na.rm = TRUE))
  }
  z[z > zlim[2]] <- zlim[2]
  levels <- pretty(zlim, 10)

  # pal <- colorRampPalette(c("lightblue", "orange"))

  if (!add){
    fields::imagePlot(x = unique.x, y = unique.y, z = z, zlim = zlim, 
                      col = faded_virdis(arg.col, min_alpha = 1),
                       asp = 1, xlim = xlim, ylim = ylim, xlab = "x", ylab = "y", legend.width = 2,
                       main = "Density surface", legend.line = 4, legend.lab = expression("Density per " ~ ha^"-1"),
                       legend.shrink = 1, legend.mar = 7)
  } else {
    image(x = unique.x, y = unique.y, z = z, zlim = zlim, 
                       col = faded_virdis(arg.col, min_alpha = 0.65), add = TRUE, legend.width = 2.5)
  }
  
  if(!is.null(session)){
    traps = get_trap(fit)[[session]]
    if(is.null(trap.plot$col)) trap_col = 1
    if(is.null(trap.plot$pch)) trap_pch = 4
    if(is.null(trap.plot$lwd)) trap_lwd = 2
    points(traps, col = trap_col, pch = trap_pch, lwd = trap_lwd)
  }
  
  if (plot.contours){
    suppressWarnings(contour(x = unique.x, y = unique.y, z = z, levels = levels,
            drawlabels = TRUE, add = TRUE, asp = 1))
  }

  invisible(pred)
}


#' Title
#'
#' @param x
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.acre_data <- function(x, ...){
  
  extra_args = list(...)
  type = extra_args$type
  session = extra_args$session
  ask = extra_args$ask
  if(is.null(ask)) ask = TRUE
  xlim = extra_args$xlim
  ylim = extra_args$ylim
  
  if(is.null(type)) stop('argument "type" is needed, which should be either "survey", "capt" or "covariates".')
  n.sessions = length(x$traps)
  if(!(type %in% c('survey', 'capt', 'covariates'))) stop('invalid input for "type", which should be either "survey", "capt" or "covariates".')
  ################################################################################################
  if(type == 'survey'){
    if(is.null(session)){
      session = 1
    }
    
    stopifnot(session <= n.sessions)
    
    masks = get_mask_from_data(x)[[session]]
    
    if(is.null(extra_args$pch_mask)) pch_mask = "." else pch_mask = extra_args$pch_mask
    if(is.null(extra_args$pch_trap)) pch_trap = 16 else pch_trap = extra_args$pch_trap
    if(is.null(extra_args$col_mask)) col_mask = 8 else col_mask = extra_args$col_mask
    if(is.null(extra_args$col_trap)) col_trap = "red" else col_trap = extra_args$col_trap
    if(is.null(extra_args$cex_mask)) cex_mask = 2 else cex_mask = extra_args$cex_mask
    if(is.null(extra_args$cex_trap)) cex_trap = 1 else cex_trap = extra_args$cex_trap
    
    plot(masks, pch = pch_mask, cex = cex_mask, asp = 1, col = col_mask)
    points(get_trap_from_data(x)[[session]], pch = pch_trap, col = col_trap, cex = cex_trap)
  }
  ################################################################################################
  if(type == 'capt'){
    # When session is not assigned, plot all capture history by default
    if(is.null(session)) session = 0
    stopifnot(session <= n.sessions)
    
    # In this "type", argument cex (cex_det) is used to extra_args the traps' symbol
    # and the cex of the circle used to show the activated traps (cex_capt) should be always cex_det + 2
    # in order to cover the symbol of detectors
    if(is.null(extra_args$cex)){
      cex_det = 3
      cex_capt = 5
    } else {
      cex_det = extra_args$cex
      if(cex_det > 7){
        warning("The size of symbol of traps might be too big to show other information properly.")
      }
      cex_capt = extra_args$cex + 2
    } 
    
    # Get capture history from data
    # Be careful that, in this function, the animal_ID and ID will be
    # converted to natural successive numbers
    capt = get_capt_for_plot(x)
    c_names = colnames(capt)
    if("bearing" %in% c_names) is.bearing = TRUE else is.bearing = FALSE
    if("dist" %in% c_names) is.dist = TRUE else is.dist = FALSE
    if("ss" %in% c_names) is.ss = TRUE else is.ss = FALSE
    if("toa" %in% c_names) is.toa = TRUE else is.toa = FALSE
    if("animal_ID" %in% c_names) animal.model = TRUE else animal.model = FALSE

    # When animal ID or call ID equals to zero, it means all detection
    if(any(!is.null(extra_args$animal_id), !is.null(extra_args$call_id)) & session == 0){
      session = 1
    }
    
    if(animal.model){
      if(is.null(extra_args$animal_id)){
        a_id = 0
        # Cannot plot with call ID but without animal_ID
        if(!is.null(extra_args$call_id)) stop("Please provide information about animal_ID.")
      } else {
        a_id = extra_args$animal_id
        stopifnot(length(a_id) == 1)
      }
    }
    
    if(is.null(extra_args$call_id)){
      c_id = 0
    } else {
      c_id = extra_args$call_id
    } 

    t_list = get_trap_from_data(x)
    m_list = get_mask_from_data(x)
    ##################################################################################
    # If session is 0, plot all sessions, create a vector for a sessions for plot
    if(session == 0) session = seq(n.sessions)
    
    for(s in session){
      traps = as.data.frame(t_list[[s]])
      masks = m_list[[s]]
      buffer = attr(masks, 'buffer')
      masks = as.data.frame(masks)
      capt_session = subset(capt, capt$session == s)
      
      if(animal.model){
        if(a_id != 0){
          capt_session = subset(capt_session, capt_session$animal_ID %in% a_id)
          if(all(c_id != 0)){
            capt_session = subset(capt_session, capt_session$ID %in% c_id)
          }
        }
        keys = paste(capt_session$animal_ID, capt_session$ID)
      } else {
        if(all(c_id != 0)){
          capt_session = subset(capt_session, capt_session$ID %in% c_id)
        }
        keys = capt_session$ID
      }
      
      if(length(keys) == 0) stop("Nothing to plot, please double check the assigned animal_ID or ID.")
      
      u_keys = unique(keys)
      
      # Because argument of x,ylim should work for each session, we cannot override "xlim" or "ylim",
      # so we create 2 tem variables xlim_plot and ylim_plot for each session
      if(is.null(xlim)){
        xlim_plot = range(masks$x)
      } else {
        xlim_plot = xlim
      }
      
      if(is.null(ylim)){
        ylim_plot = range(masks$y)
      } else {
        ylim_plot = ylim
      }
      
      # Plot calls one by one
      for(k in u_keys){
        i_k = which(keys == k)
        one_call = capt_session[i_k,,drop = FALSE]
        
        # Base plot (mask)
        plot(NA, xlim = xlim_plot, ylim = ylim_plot, asp = 1, 
             xlab = "x", ylab = "y")
        
        # Add title
        if(animal.model){
          plot_title = paste0("session: ",s, ", animal ID: ", one_call$animal_ID[1],
                              ", call ID: ", one_call$ID[1])
        } else {
          plot_title = paste0("session: ",s, ", call ID: ", one_call$ID[1])
        }
        
        title(main = plot_title)
        
        # Add traps
        points(x = traps[,1], y = traps[,2], cex = 1.2, col = "red", 
               lwd = 1.5, pch = 4)
        
        # Circle active traps
        activated_traps = traps[one_call$trap,]
        points(x = activated_traps[,1], y = activated_traps[,2], 
               cex = 2, col = "red", lwd = 1.5, pch = 1)
        
        if(is.bearing){
          if(is.dist){
            arrow_len = one_call$dist
          } else {
            if(is.null(extra_args$arrow_len)) arrow_len = 0.382 * buffer else arrow_len = extra_args$arrow_len
          }
          sinb <- sin(one_call$bearing)*arrow_len
          cosb <- cos(one_call$bearing)*arrow_len
          
          arrows.df <- data.frame(
            x = activated_traps[, 1],
            y = activated_traps[, 2],
            xend = activated_traps[, 1] + sinb,
            yend = activated_traps[, 2] + cosb
          )
          
          # Remove any NA bearing captures
          arrows.df <- arrows.df[complete.cases(arrows.df), ]
          
          arrows(arrows.df$x, arrows.df$y, arrows.df$xend, arrows.df$yend, col="red", lwd = 1.5, length = 0.07)
          # plot_one_call = plot_one_call + 
          #   geom_segment(data = activated_traps, mapping = aes(x = activated_traps$x, y = activated_traps$y,
          #                                                    xend = activated_traps$x + sin(one_call$bearing) * arrow_len,
          #                                                    yend = activated_traps$y + cos(one_call$bearing) * arrow_len),
          #                arrow = arrow(length = unit(0.02, "npc")), colour = "blue")
          
        } else if(is.dist){
          if(is.null(extra_args$circle_acc)) npoints = 50 else npoints = extra_args$circle_acc
          # Generate circle points
          data_cir = circle_fun(centre = activated_traps, r = one_call$dist, npoints = npoints)
          
          # Plot distance circle around each trap
          for (i in unique(data_cir$cir_index)) {
            cir <- subset(data_cir, cir_index == i)
            polygon(cir$x, cir$y,
                    border = "red",
                    col = NA,
                    lwd = 0.5,
                    lty = "dashed")
          }
        } 
        
        if (is.toa) {
          toa_order <- rank(one_call$toa, ties.method = "min")
          # Annotate traps
          text(activated_traps[, "x"] + 1, activated_traps[, "y"] - 1, labels = toa_order,
               col = "black", cex = 0.6, adj = c(0, 1))
        }
        
        if (ask) {
          # Make sure to only ask if we are plotting more than one plot, 
          # and it is not the last plot in the list
          if (length(u_keys) > 1 && !(k == u_keys[[length(u_keys)]] & s == session[length(session)])) {
            prompt_user_for_next_plot()
          }
        }
      }
      # End of session s 
    }
  # End of type "capt"
  }
  

  ################################################################################################

  if(type == 'covariates'){
    if(is.null(session)){
      session = 1
    }
    
    masks = as.data.frame(get_mask_from_data(x)[[session]])
    masks$mask = seq(nrow(masks))
    
    if(is.null(xlim)) xlim <- range(masks[, 'x'])
    if(is.null(ylim)) ylim <- range(masks[, 'y'])
    
    masks = subset(masks, masks$x <= xlim[2] & masks$x >= xlim[1] & masks$y <= ylim[2] & masks$y >= ylim[1])
    
    masks_mat = as.matrix(masks[, c('x', 'y'), drop = FALSE])
    D_cov_for_model = x$par.extend$data$mask

    
    if(is.null(D_cov_for_model)){
      stop('There is no location related covariates, nothing to plot.')
    }
    
    
    D_cov_for_model = D_cov_for_model[D_cov_for_model$session == session,
                                      -which(colnames(D_cov_for_model) == 'session')]
    
    if (!is.null(extra_args$select_cov)) {
      D_cov_for_model = D_cov_for_model[,which(colnames(D_cov_for_model) %in% c(extra_args$select_cov, "mask"))]
    }
    
    D_cov_for_model = merge(D_cov_for_model, masks, by = 'mask', all = FALSE)
    D_cov_for_model = D_cov_for_model[order(D_cov_for_model$mask),]
    D_cov_for_model = D_cov_for_model[, -which(colnames(D_cov_for_model) == 'mask')]
    
    cov_list = colnames(D_cov_for_model)
    cov_list = cov_list[-which(cov_list == 'x' | cov_list == 'y')]
    
    
    unique.x <- sort(unique(masks[, 'x']))
    unique.y <- sort(unique(masks[, 'y']))
    
    if(is.null(extra_args$arg.col)) arg.col = 100 else arg.col = extra_args$arg.col
    
    if(is.null(extra_args$plot.contours)){
        plot.contours = FALSE
    } else {
        plot.contours = extra_args$plot.contours
    }
    
    for(i in cov_list){
      if(is(D_cov_for_model[[i]], 'numeric')){
        
        z <- squarify(masks_mat, D_cov_for_model[[i]])
        zlim <- range(z, na.rm = TRUE)
        
        fields::imagePlot(x = unique.x, y = unique.y, z = z, zlim = zlim, col = faded_virdis(arg.col, min_alpha = 1),
                           asp = 1, xlim = xlim, ylim = ylim, xlab = "x", ylab = "y",
                           main = paste0("Plot of covariate ", i, ", for session ", session),
                           legend.width = 2, legend.shrink = 1, legend.mar = 7)
        
        if(plot.contours){
          levels <- pretty(zlim, 10)
          contour(x = unique.x, y = unique.y, z = z, levels = levels, drawlabels = TRUE, add = TRUE)
        }
      } else {
        
        v = D_cov_for_model[[i]]
        v_num = as.numeric(as.factor(v))
        
        match_table = data.frame(code = v_num, name = v)
        match_table = match_table[!duplicated(match_table),,drop = FALSE]
        
        z <- squarify(masks_mat, v_num)
        zlim <- range(match_table$code)
        
        col = faded_virdis(arg.col, min_alpha = 1)
        idx = round(scale_convert(match_table$code, seq(arg.col)),0)
        
        col_leng = col[idx]
        
        image(x = unique.x, y = unique.y, z = z, zlim = zlim, col = col, 
              asp = 1, xlim = xlim, ylim = ylim, xlab = "x", ylab = "y",
              main = paste0("Plot of covariate ", i, ", for session ", session))
        legend(x = "topright", legend = match_table$name, fill = col_leng)
        
      }
      
      if (ask) {
        # Make sure to only ask if we are plotting more than one plot, 
        # and it is not the last plot in the list
        if (length(cov_list) > 1 && i != cov_list[[length(cov_list)]]) {
          prompt_user_for_next_plot()
        }
      }
    }
    #end of type == 'covariates'
  }
  
    
}


#' Title
#'
#' @param x 
#' @param ... 
#' 
#' @return
#' @export
#'
#' @examples
plot.acre_tmb = function(x, ...){
  extra_args = list(...)
  type = extra_args$type
  
  if(is.null(type)){
    stop('argument "type" is needed, which should be either "survey", "capt", 
         "covariates", "detfn", "locations" or "Dsurf".')
  }
  
  if(type %in% c('survey', 'capt', 'covariates')){
    plot.acre_data(x = x$args, ...)
  } else if(type == 'detfn'){
    args_pass = list(fit = x, ...)
    args_pass$type = NULL
    do.call('show.detfn', args_pass)
    
  } else if(type == 'Dsurf'){
    args_pass = list(fit = x, ...)
    args_pass$type = NULL
    do.call('show.Dsurf', args_pass)

  } else if(type == 'locations'){
    args_pass = list(fit = x, ...)
    args_pass$type = NULL
    do.call('plot_locations', args_pass)
    
  } else if(type == 'dev') {
    args_pass = list(fit = x, ...)
    args_pass$type = NULL
    do.call('plot_dev', args_pass)
    
  } else {
    stop('invalid input for "type", which should be either "survey", "capt", 
         "covariates", "detfn", "locations" or "Dsurf".')
  }
}


prompt_user_for_next_plot <- function() {
  if (isRStudio <- Sys.getenv("RSTUDIO") == "1") {
    input <- readline(prompt = "Hit <Return> to see next plot or <Esc> to cancel:")
  } else {
    input <- readline(prompt = "Hit <Return> to see next plot or <Ctrl-c> to cancel:")
  }
}

# Empty plot with mask limits and traps
# Call with traps = NULL explicitly to plot without traps
plot_dev <- function(fit, session = 1, 
                     mask = get_mask(fit)[[session]], 
                     traps = get_trap(fit)[[session]],
                     call_id = NULL, ...) {
  plot(NA, 
       xlim = range(mask[,1]), 
       ylim = range(mask[,1]),
       xlab = "x", ylab = "y",
       ...
       )
  grid(col = "grey85", lty = 1, lwd = 0.7)
  
  if (!is.null(traps)) {
    plot_traps(fit, call_id = call_id)
  }
}

plot_traps <- function(fit, call_id=NULL, animal_id=NULL, session=1) {
  traps <- get_trap(fit)[[session]]
  
  # Plot the trap array
  points(x = traps[,1], y = traps[,2], cex = 1.2, col = "red", 
         lwd = 1.5, pch = 4)
  
  if (!is.null(call_id)) {
    if (length(call_id) > 1 | !is.numeric(call_id)) {
      stop(paste("'call_id' must be a single numeric value.", 
                 "Plotting traps is only supported for single calls."))
    }
    
    # Grab appropriate capture data
    bincapt <- get_capt_by_id(fit, call_id, animal_id, session, 
                              return_bincapt = T)[,-1]
    
    active_traps <- subset(traps, bincapt == 1)
    # Circle the activated traps
    points(x = active_traps[,1], y = active_traps[,2], 
           cex = 2, col = "red", lwd = 1.5, pch = 1)
  }
}

plot_bearing_arrows <- function(fit, call_id, animal_id=NULL, session=1) {
  # Make sure call_id is a single numeric
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting bearing arrows is only supported for single calls."))
  }
  
  if (!("bearing" %in% fit$infotypes)) {
    stop("'fit' does not contain bearing information.")
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)[,-1]
  
  is.dist <- "dist" %in% fit$infotypes
  traps <- get_trap(fit)[[session]]
  activated_traps <- subset(traps, bincapt == 1)
  bearing <- capt$bearing
  dist <- capt$dist
  
  # Set appropriate arrow length. Default to 1/3 buffer radius.
  arrow_len <- rep(0.333 * get_buffer(fit)[session], length(capt$bearing))
  if (is.dist) {
    # Account for the possibility that dist may be NA for some captures
    arrow_len[!is.na(capt$dist)] <- capt$dist[!is.na(capt$dist)]
  }
  
  sinb <- sin(bearing)*arrow_len
  cosb <- cos(bearing)*arrow_len
  
  arrows.df <- data.frame(
    x = activated_traps[, 1],
    y = activated_traps[, 2],
    xend = activated_traps[, 1] + sinb,
    yend = activated_traps[, 2] + cosb
  )
  
  # Remove any NA bearing captures
  arrows.df <- arrows.df[complete.cases(arrows.df), ]
  
  arrows(arrows.df$x, arrows.df$y, arrows.df$xend, arrows.df$yend, col="red", lwd = 1.5, length = 0.07)
}

plot_distance_circles <- function(fit, call_id, animal_id=NULL, session=1) {
  # Make sure call_id is a single numeric
  if (length(call_id) > 1 | !is.numeric(call_id)) {
    stop(paste("'call_id' must be a single numeric value.", 
               "Plotting distance circles is only supported for single calls."))
  }
  
  if (!("dist" %in% fit$infotypes)) {
    stop("'fit' does not contain distance information.")
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  bincapt <- get_capt_by_id(fit, call_id, animal_id, session, return_bincapt=T)[,-1]
  dist <- capt$dist
  traps <- get_trap(fit)[[session]]
  activated_traps <- subset(traps, bincapt == 1)
  
  # Remove NA distances
  keep <- !is.na(capt$dist)
  # Note the drop=FALSE keeps the data structure; makes easier to loop over
  activated_traps <- activated_traps[keep, ,drop = FALSE]
  dists <- capt$dist[keep]
  
  for (i in seq_len(nrow(activated_traps))) {
    
    coords <- circle_fun(centre = activated_traps[i, ,drop = FALSE],
                         r = dists[i])
    
    polygon(coords$x, coords$y,
            border = "red",
            col = NA,
            lwd = 0.5,
            lty = "dashed")
  }
}

plot_toa_order <- function(fit, call_id, animal_id=NULL, session=1) {
  if (length(call_id) != 1L || !is.numeric(call_id)) {
    stop("'call_id' must be a single numeric value; ",
         "plotting TOA order is only supported for a single call.")
  }
  
  if (!("toa" %in% fit$infotypes)) {
    stop("'fit' does not contain time-of-arrival information.")
  }
  
  # Grab appropriate capture data
  capt <- get_capt_by_id(fit, call_id, animal_id, session)
  
  # Remove NA distances
  capt <- capt[complete.cases(capt), ]
  
  toa <- capt$toa
  
  # Figure out the order of arrival
  traps <- get_trap(fit)[[session]][capt$trap, ,drop = FALSE]
  toa_order <- rank(toa, ties.method = "min")

  # Annotate traps
  text(traps[, 1] + 1, traps[, 2] - 1, labels = toa_order,
       col = "black", cex = 0.6, adj = c(0, 1))
}



