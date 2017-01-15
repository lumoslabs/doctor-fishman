
## multi-arm split test analysis
options(shiny.maxRequestSize = 100 * 1024^2)

#####################
# DATA UPLOAD FUNCTIONS
#####################

# upload rate data
load_test_data <- function(f) {
  require(data.table)
  if(is.null(f)) return(NULL)
  else {
    dt <- data.table(read.csv(f, na.strings = c("NA", "UNKNOWN", "NotApplicableToUser")))
    # Note, downstream in calculate_aggregates(), an error will occur if the trials and successes column names are
    # equivalent to the name of the variable in that function (stupid problem)
    # Here we greedily check for those "reservered" words and adjust them so that they will not cause an error
    # downstream.
    if("trials_label" %in% names(dt)) setnames(dt, "trials_label", "trials__label")
    if("successes_label" %in% names(dt)) setnames(dt, "successes_label", "successes__label")
    return(dt)
  }
}

typed_in_data_is_formatted <- function(num_groups, input) {
  require(stringr)
  # before we attempt to extract the typed-in data, lets make sure that it is in the expected format
  two_commas = all(sapply(1:num_groups, function(idx) {
    mytext = input[[paste('g',idx,'_data', sep = "")]]
    str_count(mytext, ',') == 2 # we expect 2 commas for each input text
  }))

  # TODO: write additional criteria
  #       verify that second input is number and third input is number
  #       verify that there are no blanks / NAs

  num_criteria = 2
  criteria = c(two_commas, TRUE)

  all((length(criteria) == num_criteria), # no NULL crept in?
      !is.na(criteria),  #no NA crept in?
      criteria) # No outright falses?

}

extract_typed_data <- function(num_groups, input) {
  require(data.table)
  output = rbindlist(lapply(1:num_groups, function(idx) {
    mytext = input[[paste('g',idx,'_data', sep = "")]]
    read.csv(text = mytext, header = FALSE, col.names = c("arm", "trials", "successes"))
  }))
  output[, rate := (successes / trials) * 100]
  output
}

#####################
# CONVERSION RATE FUNCTIONS
#####################

# calculates aggregate values
calculate_aggregates <- function(dt, split, trials_label, successes_label) {
  if(split %in% names(dt)) setnames(dt, split, "arm")

  aggs <- dt[, list(trials = sum(get(trials_label)), successes = sum(get(successes_label))), by = arm]

  aggs <- aggs[arm != "NotApplicableToUser" & !is.na(arm)]
  aggs[, rate := (successes / trials) * 100]
  aggs <- aggs[, list(arm, trials, successes, rate)]
  return(aggs)
}

# draw samples from beta distribution
get_samples <- function(params, n = 1000000) {
  samples <- mapply(function(alpha, beta) rbeta(n, alpha, beta), params$alpha, params$beta)
  samples = 100 * samples # % instead of fraction
  return(samples)
}

# find the losing test arm for each sample
get_losers <- function(samples) {
  loser <- max.col(-samples)
  return(loser)
}

# calculate the proportion of samples where each arm wins
get_ps <- function(num_arms, vals) {
  freqs <- data.table(arm = vals)[, list(n = .N), by = "arm"]
  freqs <- merge(data.table(arm = 1:num_arms), freqs, by = "arm", all.x = T)
  setkey(freqs, arm)
  freqs[is.na(n), n := 0]
  return(freqs[, n / sum(n)])
}

# build an analysis table based on sampling from the posterior beta distributions
analyze <- function(dt, samples_dt, progress_obj=NULL) {
  require(plyr)
  if(is.null(dt)) return(NULL)
  else {
    params <- data.table(arm = dt$arm, alpha = dt$successes + 1, beta = dt$trials - dt$successes + 1)
    if(!is.null(progress_obj)) progress_obj$inc(1/5, "sampling from beta")
    samples <- get_samples(params)
    if(!is.null(progress_obj)) progress_obj$inc(1/5, "calculating winners and losers for samples")

    # Get Prob of Worst Loser
    losers <- get_losers(samples)
    if(!is.null(progress_obj)) progress_obj$inc(1/5, "calculating win/lose probabilities")
    p_worsts <- get_ps(length(dt$arm), losers)
    if(!is.null(progress_obj)) progress_obj$inc(1/5, "calculating confidence intervals")
    CI_lower <- apply(samples, 2, function(x) quantile(x, .025))
    CI_upper <- apply(samples, 2, function(x) quantile(x, .975))

    # this table will only report simulations where there was a solo winner (ties are reported separately)
    if(!is.null(progress_obj)) progress_obj$inc(1/5, "calculating win XOR tie probabilities")
    percents = samples_dt[tie==FALSE, list(prob_win = 100*(.N)/nrow(samples_dt)), by = c("winner_idx")]

    if(!is.null(progress_obj)) progress_obj$inc(1/5, "Formatting results table")
    variant_list = as.data.table(data.frame(Variant = dt$arm, winner_idx = c(1:length(dt$arm))))
    percents = merge(percents, variant_list, by = "winner_idx", all.y = TRUE)
    percents[is.na(prob_win), prob_win:= 0]

    result <- data.table(arm = params$arm,
                         `Lower 95% CI` = CI_lower,
                         `Upper 95% CI` = CI_upper,
                         `Probability of worst loser (%)` = p_worsts * 100)

    # format table for dashboard
    formatted_table = copy(dt)
    setnames(formatted_table, "successes", "Successes")
    setnames(formatted_table, "trials", "# Observations")
    setnames(formatted_table, "rate", "rate (%)")
    formatted_table = merge(formatted_table, result, by = c("arm"))
    setnames(formatted_table, "arm", "Variant")
    formatted_table = merge(percents, formatted_table, by = c("Variant"))
    setnames(formatted_table, "prob_win", "Probability of winner (%)")
    desired_col_order = c("Variant", "Probability of winner (%)", "# Observations", "Successes", "rate (%)",
                          "Lower 95% CI", "Upper 95% CI", "Probability of worst loser (%)")
    formatted_table = formatted_table[, desired_col_order, with = F]
    return(formatted_table)
  }
}

diff_function <- function(variant, baseline){
  diff = (variant - baseline) / baseline * 100.0
  diff[diff==Inf] = 1000 # edge case handling
  return(diff)
}

calculate_differences <- function(dt, baseline_variant, do_progress_bar = TRUE) {
  require(data.table)
  if(do_progress_bar){
    progress <- shiny::Progress$new()
    progress$inc(1/5, "calculating difference distributions")
    on.exit(progress$close())
  }
  if(tolower(baseline_variant) %in% tolower(dt$arm)){
    params <- data.table(arm = dt$arm, alpha = dt$successes + 1, beta = dt$trials - dt$successes + 1)
    samples <- as.data.table(get_samples(params, 100000))
    setnames(samples, as.character(dt$arm))
    diff_dt = as.data.table(lapply(samples, function(x) diff_function(x, samples[[baseline_variant]])))
    return(diff_dt)
  } else return(NULL)
}

find_which_simulations_were_ties <- function(samples_dt, this_winner, min_effect) {
  # This function adds a column to samples_dt which is logical for whether or not this row represents a "tie" for
  # winner, given the min_effect.
  # Note: this is a dummy function created for a specific apply statement

  # We expect samples_dt to be a data.table which has the following columns:
  #     1...n :: where n is the number of variants in the test, these can be named whatever
  #     sim_idx :: refers to which simulation this is from the entire set of samples.
  #     winner_idx :: refers to which column 1...n has the largest value for each row.
  # samples_dt is a subset of the rows returned from get_samples(), in which one specific variant is the winner in all
  # the rows, i.e., length(unique(samples_dt$winner_idx)) = 1

  # get subset of simulations where this variant won
  subset_samples = samples_dt[winner_idx == this_winner]
  winning_variant = names(subset_samples)[this_winner]
  setnames(subset_samples, winning_variant, "winner_")

  # determine the "tie" range (i.e., win_criteria) for all of these simulations
  subset_samples[, win_criteria:= winner_ - winner_ * min_effect / 100]
  candidates = subset_samples[, !names(subset_samples) %in% c("winner_"), with = F]

  # find next best winner and see if it is within the win_criteria
  next_best = as.data.frame(candidates[, c(!names(candidates) %in% c(winning_variant, "sim_idx", "win_criteria", "winner_idx")), with = F])
  next_best_idx <-max.col(next_best, ties.method='first')  #
  candidates[, next_best_value:= next_best[cbind(1:nrow(next_best), next_best_idx)]]
  candidates[, tie:= next_best_value > win_criteria ]

  return(merge(samples_dt, candidates[, c("tie", "sim_idx", "win_criteria"), with = F], by = c("sim_idx")))
}

make_recommendation_message <- function(analysis_table, danger_threshold, min_effect_size, tied_probability, aggs) {
  # This applies a decision criteria and generates a plain language recommendation of actions to take given the data submitted to the dashboard

  # parameters of the decision criteria
  value_remaining_criteria = 2 # if the worst case scenario is a loss by 2% or less, then we call it a winner
  losering_criteria = 10 # if something has a probabilty of winning < losering_criteria, we recommend ditching it

  true_winner = as.character(analysis_table[which.max(analysis_table[["rate (%)"]])]$Variant)
  true_value_remaining = get_value_remaining(aggs, true_winner, danger_threshold)
  #given true_value_remaining, tied_probability, value_remaining_criteria, min_effect_size, losering_criteria
  if(true_value_remaining<=(max(value_remaining_criteria, min_effect_size))){
    # We are ready to winner
    line2_1 = p("Given your data, your minimum effect size of interest, and your worst-case scenario tolerance, you should:")
    line2_2 = p(strong(paste("Call the winner:", true_winner)))
    recommendation_message = paste(line2_1, line2_2, sep = '<br/>')
  } else if(any(tied_probability[["prob_win"]]<losering_criteria)){
    # Lose some and keep collecting
    line2_1 = p("Given your data, your minimum effect size of interest, and your worst-case scenario tolerance, you should:")
    line2_2 = p(strong("1)"), " LOSER the following variants: ", paste(as.character(tied_probability[tied_probability[["prob_win"]]<losering_criteria]$Variant), collapse = ", "))
    line2_3 = p(strong("2)"), "Collect more data for the remaining variants.")
    recommendation_message = paste(line2_1, line2_2, line2_3, sep = '<br/>')
  } else{
    # collect more data
    line2_1 = p("Given your data, your minimum effect size of interest, and your danger tolerance, you should:")
    line2_2 = p(strong("Collect more data."))
    recommendation_message = paste(line2_1, line2_2, sep = '<br/>')
  }
  return(recommendation_message)

}


#####################
# BOOTSTRAPPING FUNCTIONS
#####################

get_boot_results <- function(dt, groupby_col, metric_col, progress_obj) {
  require(parallel)
  groups = unique(as.character(dt[[groupby_col]]))
  data_list = lapply(groups, function(x) dt[(dt[[groupby_col]]== x)][[metric_col]])
  names(data_list) = groups
  progress_obj$inc(1/5, paste("bootstrapping for ", groups, " in parallel.", sep = ""))
  boot_results = mclapply(data_list, do_boot_wrapper, mc.cores = min(5,length(groups))) # theoretically could set mc.cores = min(12,length(groups)) here for higher performance, but risks the possibility of 144 processes
  return(boot_results)
}

do_boot_wrapper <- function(data) {
  boot_iter = 10000
  b = mclapply(1:boot_iter, function(x) {
    mean(data[sample(1:length(data), size = length(data), replace=T)])
  }, mc.cores = 12)
  return(unlist(b))
}

extract_ci_from_boot_results <- function(b, variant_name) {
  ci_percent = 0.95
  lconf = quantile(b, c((1 - ci_percent) / 2))
  uconf = quantile(b, c(1 - (1 - ci_percent) / 2))
  return(data.table(variant=variant_name, l.conf=lconf, u.conf=uconf))
}

get_prob_of_winner_from_boot_results <- function(boot_results){
  variant_list = as.data.table(data.frame(variant = names(boot_results),
                                     winner_idx = c(1:length(names(boot_results)))
                                     )
                          )
  boot_results[, winner_idx:= max.col(boot_results)]
  percents = boot_results[, list(prob_win = 100*(.N)/nrow(boot_results)), by = c("winner_idx")]
  percents = merge(percents, variant_list, by = "winner_idx", all.y = TRUE)
  percents[is.na(prob_win), prob_win:= 0]
  percents[, winner_idx:=NULL]
  return(percents)
}

run_bootstrapping_steps <- function(numeric_data, progress_obj){
  boot_results = get_boot_results(numeric_data, "variant", "data_value", progress_obj)
  progress_obj$inc(1/3, "extracting confidence intervals")
  boot_summary = rbindlist(lapply(c(1:length(boot_results)), function(x) extract_ci_from_boot_results(boot_results[[x]], names(boot_results)[[x]]) )  ) # this is a verbose way to extract the specific confidence interval info from each boot obj into a data.table
  progress_obj$inc(1/3, "getting prob(winner)")
  prob_winner_stats = get_prob_of_winner_from_boot_results(as.data.table(boot_results))

  # format the bootstrap output table
  data_agg = numeric_data[, list(counts = .N,
                              mean = mean(data_value, na.rm = T)),
                       by = c("variant")]
  data_summary = merge(data_agg, boot_summary, by = c("variant"))
  data_summary = merge(data_summary, prob_winner_stats, by = c("variant"))

  #re order for consistency across plots
  data_summary$variant = as.character(data_summary$variant)
  setorderv(data_summary, "variant") # order by alphabetical of variant name for consistency across all plots
  return(data_summary)
}

format_data_summary <- function(data_summary){
  # format table
  setnames(data_summary, "variant", "Variant")
  setnames(data_summary, "prob_win", "Probability of winner (%)")
  setnames(data_summary, "counts", "# Observations")
  setnames(data_summary, "l.conf", "Lower 95% CI")
  setnames(data_summary, "u.conf", "Upper 95% CI")

  output_order = c("Variant", "Probability of winner (%)", "# Observations", "mean", "Lower 95% CI", "Upper 95% CI")
  data_summary = data_summary[, output_order, with = F]
  return(data_summary)
}

#####################
# VALUE REMAINING FUNCTIONS
# calculations taken from: https://support.google.com/analytics/answer/2846882?vid=1-635778645345403433-155365758
#####################

get_value_remaining <- function(aggs, candidate_winner, danger_scenario_threshold){
  # returns 95%tile of value remaining, as a percent (%)
  # utilizes the pre-existing calculate_differences() function
    if(!is.null(candidate_winner) && candidate_winner != ""){
    diff_dt = calculate_differences(aggs, candidate_winner, do_progress_bar = F)

    # compute value remaining
    diff_df = as.data.frame(diff_dt)
    indx <- max.col(diff_df)
    value_remaining = diff_df[cbind(1:nrow(diff_df), indx)]
    percentile = 1 - danger_scenario_threshold / 100
    value_remaining_stat = round(quantile(value_remaining, percentile), 2)
    return(value_remaining_stat)
  } else NULL
}

#####################
# PLOTTING FUNCTIONS
#####################

# construct a plot with samples from the rate data
make_plot <- function(dt) {
  require(reshape2)
  require(data.table)
  require(ggplot2)
  require(plyr)
  if(is.null(dt)) return(NULL)
  else {
    theme_set(theme_bw(base_size=18))
    params <- data.table(arm = dt$arm, alpha = dt$successes + 1, beta = dt$trials - dt$successes + 1)
    samples <- data.table(get_samples(params, 100000))
    setnames(samples, names(samples), as.character(dt$arm))

    # this ensures that they are alphabetical, so that order will match other plots.
    setcolorder(samples, sort(names(samples)))
    samp_dt <- melt(data = samples, measure.var = names(samples), variable.name = "arm", value.name = "p")
    p <- qplot(p, data = samp_dt, color = arm, geom = "density", size = I(2)) +
      xlab("rate (100*successes/trials)") + ylab("density")
    return(p)
  }
}

# differences in rate plot
make_differences_plot <- function(dt, baseline_variant) {
  require(data.table)
  require(ggplot2)
  if(is.null(dt)) return(NULL)
  else {
    if(tolower(baseline_variant) %in% tolower(dt$arm)){
      diff_dt = calculate_differences(dt, baseline_variant)
      diff_data_for_plot = rbindlist(lapply(c(1:length(diff_dt)),
                                      function(x) extract_ci_from_boot_results(diff_dt[[x]], names(diff_dt)[[x]]) )  )
      diff_data_for_plot[, mean:= as.numeric(lapply(c(1:length(diff_dt)), function(x) mean(diff_dt[[x]], na.rm =T)))]
      return(make_nice_bar_plot(diff_data_for_plot, group_var = "variant", y = "mean", up = "u.conf",
                                lw = "l.conf", ylab = "% Difference",
                                plot_title = paste("% Delta from", baseline_variant,"with 95% CI")))
    } else{
      return(NULL)
    }
  }
}

truncate_labels <- function(str_vector, truncate_limit = 10){
  # This takes a vector of character strings and returns vector of character strings with truncated titles,
  # for the purposes of plotting and displaying labels.
  trunc_labels = substring(str_vector, 1, truncate_limit)
  index = nchar(str_vector)>truncate_limit
  append_vector = rep("", length(str_vector))
  append_vector[index] = "..."
  trunc_labels = paste(trunc_labels, append_vector, sep = "")
  return(trunc_labels)
}

make_nice_bar_plot <- function(summary_data, group_var = 'variable', y = 'my_mean', up = "up", lw = "lw", xlab = "",
                               ylab = "", ylim_ = NULL, plot_title = "", legend_position = "right") {
  require(data.table)
  require(ggplot2)
  if(is.null(summary_data)) return(NULL)
  else {
    # get plot limits
    if(is.null(ylim_)){
      yrange = max(summary_data[[up]]) - min(summary_data[[lw]])
      ylim_min = min(summary_data[[lw]]) - 0.1 * yrange
      ylim_max = max(summary_data[[up]]) + 0.1 * yrange
    } else {
      ylim_min = ylim_[1]
      ylim_max = ylim_[2]
    }

    p = ggplot(summary_data, aes_string(x = group_var, y = y, fill = group_var)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_x_discrete(labels=sort(truncate_labels(summary_data[[group_var]]))) +
      geom_errorbar(aes_string(ymax = up, ymin = lw),
                    position = position_dodge(0.5),
                    data = summary_data,
                    colour = "black") +
      theme_bw() +  coord_cartesian(ylim=c(ylim_min, ylim_max)) +
      theme(axis.text.x = element_text(size = 17, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
            axis.text.y = element_text(size = 17, hjust = 0.5, color = "black"),
            axis.title.y = element_text(size = 17, hjust = 0.5),
            plot.title = element_text(size = 17, vjust = 2)) +
      xlab(xlab) + ylab(ylab) + ggtitle(plot_title)  +
      theme(legend.position = legend_position)
    return(p)
  }
}

#####################
# SHINY SERVER FUNCTION
#####################

shinyServer(function(input, output, session) {

  # this will monitor the button clicks from either upload or type-in mode.
  myReactVals <- reactiveValues(run_data_button = 0,
                                doPlot = FALSE,
                                danger_criteria = 5,
                                baseline_variant = "control",
                                candidate_winner = NULL)

  observeEvent(input$run_data_button, {
    myReactVals$doPlot = TRUE
    myReactVals$run_data_button = myReactVals$run_data_button + 1
  })

  observeEvent(input$data_input_type, {
    # the plot disappears, and defaults reset if you flip to a different input type
    myReactVals$doPlot = FALSE
    myReactVals <- reactiveValues(run_data_button = 0,
                                  doPlot = FALSE,
                                  danger_criteria = 5,
                                  baseline_variant = "control",
                                  candidate_winner = NULL)
    updateSelectInput(session, "baseline_variant1_input", choices = "control", selected="control")
    updateSelectInput(session, "candidate_winner", choices = "", selected = NULL)
    updateNumericInput(session, "danger_criteria1", value = 5)
    updateNumericInput(session, "min_effect_size1", value = 0)
  })

  observeEvent(input$danger_criteria1, {
    myReactVals$danger_criteria = input$danger_criteria1
  })

  observeEvent(input$min_effect_size1, {
    myReactVals$min_effect_size = input$min_effect_size1
  })

  observeEvent(input$baseline_variant1_input, {
    if(is.null(input$baseline_variant1_input)){
      myReactVals$baseline_variant = "control"
    } else {
      myReactVals$baseline_variant = input$baseline_variant1_input
    }
  })

  observeEvent(input$candidate_winner, {
    myReactVals$candidate_winner = input$candidate_winner
  })

  # dynamic input fields for the type-in-data mode of the dashboard
  output$dynamic_input_fields <- renderUI({
      if(input$num_groups<=0){
        empty = "" # this shows nothing until input$num_groups >0
      } else {
        output_to_render = lapply(1:input$num_groups, function(idx) {
          textInput(paste('g', idx, '_data', sep = ""),
                    label = h5(paste("Group ", idx, " enter:  name, trials, success", sep = "")), value = NA)
        })
      }
  })

  testfilepath = reactive({
    # Define which upload path the dashboard should listen to TODO: re-factor so that we don't need redundantly name
    # elements on the ui.R side (?)
    if(input$data_input_type == 'upload') path = "testfile"
    else if(input$data_input_type == 'typein') path = "testfile"
    else if(input$data_input_type == 'upload_cash') path = "testfile_cash"
    else if(input$data_input_type == 'upload_survey') path = "testfile_survey"
    return(path)
  })

  observe({
    if(!is.null(input[[testfilepath()]])) {
      # update the drop down menu with the column names
      col_names <- names(load_test_data(input[[testfilepath()]]$datapath))
      if(input$data_input_type == 'upload' || input$data_input_type == 'typein'){
        updateSelectInput(session, "split_col", choices = col_names, selected=NULL)
        updateSelectInput(session, "trials_col", choices = col_names, selected=NULL)
        updateSelectInput(session, "successes_col", choices = col_names, selected=NULL)
      } else if(input$data_input_type == 'upload_cash') {
        updateSelectInput(session, "split_col_cash", choices = col_names, selected=NULL)
        updateSelectInput(session, "user_id_col", choices = col_names, selected=NULL)
        updateSelectInput(session, "cash_col", choices = col_names, selected=NULL)
      }else if(input$data_input_type == 'upload_survey') {
        updateSelectInput(session, "split_col_survey", choices = col_names, selected=NULL)
        updateSelectInput(session, "user_id_col_survey", choices = col_names, selected=NULL)
        updateSelectInput(session, "survey_col", choices = col_names, selected=NULL)
      }
    }
  })

  get_continuous_var_data = reactive({
    progress = shiny::Progress$new()
    on.exit(progress$close())
    progress$inc(1/3, "re-loading data")
    numerical_data = load_test_data(input[[testfilepath()]]$datapath)
    if(is.null(numerical_data)) return(NULL)
    else {
      #TODO check for users who saw multiple variants and exclude them (and report it to user)
      #TODO check for missing/NULL data explicitly, exclude it and report it to user
      #data prep
      if(input$data_input_type == "upload_cash"){
        setnames(numerical_data, input$cash_col, "data_value")
        setnames(numerical_data, input$split_col_cash, "variant")
        setnames(numerical_data, input$user_id_col, "user_id")
        numerical_data[is.na(data_value), data_value:= 0] # nulls in the data_value column are treated as non converted, so we assign them $0 cash collected
        numerical_data[, converted:= (data_value > 0)]
      } else if(input$data_input_type == "upload_survey"){
        setnames(numerical_data, input$survey_col, "data_value")
        setnames(numerical_data, input$split_col_survey, "variant")
        setnames(numerical_data, input$user_id_col_survey, "user_id")
        numerical_data = numerical_data[!is.na(data_value)]
      }
      return(numerical_data)
    }
  })

  #####################
  # RATE ANALYSIS FUNCTIONS
  #####################

  aggregates <- eventReactive(myReactVals$run_data_button, {
    if(input$data_input_type == 'upload'){
      if(!is.null(input[[testfilepath()]])) {
        output$instructions <- renderText(NULL)
        calculate_aggregates(load_test_data(input[[testfilepath()]]$datapath), input$split_col, input$trials_col, input$successes_col)
      } else {
        NULL
      }
    } else if(input$data_input_type == 'typein'){
      if(typed_in_data_is_formatted(input$num_groups, input)){
        extract_typed_data(input$num_groups, input)
      } else {
        output$instructions <- renderText("Please type in your data and format correctly")
        NULL
      }
    } else if(input$data_input_type == 'upload_cash') {
      dt = copy(get_continuous_var_data())
      dt[, trials:= 1]
      calculate_aggregates(dt, "variant", "trials", "converted")
    } else if(input$data_input_type == 'upload_survey') {
      return(NULL)
    }
  })

  aggregates_wrapper <- reactive({
    if(myReactVals$doPlot){
      aggregates()
    } else {
      NULL
    }
  })

  analysis_table <- reactive({
    progress <- shiny::Progress$new()
    progress$set(message = 'analyzing results', value = 0)
    on.exit(progress$close())
    analyze(aggregates_wrapper(), get_samples_with_ties(), progress_obj = progress)
  })

  tied_table <- reactive({
    require(reshape2)
    # this takes the simulation data from get_samples_with_ties() and aggregates them into summary info
    # TODO: This is the slow step of the dashboard. Are there ways to make this more efficient/faster?
    samples_dt = get_samples_with_ties()
    if(!is.null(samples_dt)){

      progress_obj <- shiny::Progress$new()
      progress_obj$set(message = 'computing ties', value = 0)
      on.exit(progress_obj$close())
      variant_names = as.character(aggregates_wrapper()$arm)

      #### Get specific tie results probabilities
      tie_results = samples_dt[tie==TRUE]
      if(nrow(tie_results)>0){

        tie_results[, num_winners_tied:= 0]
        for(variant in variant_names){
          tie_results[["num_winners_tied"]] = tie_results$num_winners_tied
                                                + as.numeric(tie_results[[variant]]>tie_results$win_criteria)
        }

        tie_results_melted = melt(tie_results, measure.vars = c(variant_names),
                                  id.vars = c("sim_idx", "num_winners_tied"), stringsAsFactors = FALSE)
        tie_results_melted$variable = as.character(tie_results_melted$variable)
        progress_obj$inc(1/5, "melting and ordering tie results")
        setorderv(tie_results_melted, c("sim_idx", "value"))
        num_variants = length(variant_names)
        progress_obj$inc(1/5, "identifying all tie outcomes (this may take up to 30 seconds)")

        # TODO this is the slow step because we group by the sim_idx.
        # I could parallelize it? Maybe there is a different approach?
        # This is the kind of line that makes PHP look like a great alternative - Marc
        tie_groups = tie_results_melted[, list(group=paste(sort(variable[(num_variants+1-num_winners_tied[[1]]):num_variants]), collapse = ',')), by = c("sim_idx")] # this aggregates every simulation by the tied winners

        tie_percents = tie_groups[, list(prob_win = 100*(.N)/nrow(samples_dt)), by = c("group")]
        setnames(tie_percents, "group", "Tied Winning Variants")

        inclusion_criteria = 1 # (%) a tie group must occur with this frequency to be listed in the dashboard.
        tie_percents = tie_percents[prob_win>inclusion_criteria] # filter out small contributors
        setorderv(tie_percents, "prob_win", -1)
        if(nrow(tie_percents)>0){
          setnames(tie_percents, "prob_win", paste("% chance that these variants win and are within",
                   myReactVals$min_effect_size, "% of eachother"))
          return(tie_percents)
        } else NULL
      } else NULL
    } else NULL
  })

  get_samples_with_ties <- reactive({
    # this function calls get_samples() and adds 2 additional columns:
    #    1) "tie" (logical), is the winner and 2nd place winner within the range of a tie (defined by the user-set
    #       min_effect_size)
    #    2) "win_criteria" the range within the winner, given the min_effect_size, for the next best variant to be
    #       considered a tie.
    aggs = aggregates_wrapper()
    if(!is.null(aggs)){

      # posterior for each variant
      params <- data.table(arm = aggs$arm, alpha = aggs$successes + 1, beta = aggs$trials - aggs$successes + 1)
      samples_dt = as.data.table(get_samples(params))
      variant_names = as.character(aggs$arm)

      setnames(samples_dt, variant_names)
      samples_dt[, winner_idx:= max.col(samples_dt)]
      samples_dt[, sim_idx:= rep(1:nrow(samples_dt))] # this just helps with debugging

      min_effect = 1.0 * myReactVals$min_effect_size

      # we identify ties, by subsetting the total simulations into groups of each winner.
      winner_list = unique(as.numeric(samples_dt$winner_idx))
      results = rbindlist(lapply(winner_list,
                          function(x) find_which_simulations_were_ties(copy(samples_dt), x, min_effect)))

      return(results)
    } else return(NULL)
  })

  observe({
    at = analysis_table()
    if(!is.null(analysis_table())){
      updateSelectInput(session, "candidate_winner", choices=as.character(unique(at$Variant)),
                        selected = as.character(at[which.max(at[["rate (%)"]])]$Variant))
      if("control" %in% as.character(unique(at$Variant))) {
        updateSelectInput(session, "baseline_variant1_input", choices=as.character(unique(at$Variant)), selected="control")
      } else {
        updateSelectInput(session, "baseline_variant1_input", choices=as.character(unique(at$Variant)), selected=NULL)
      }
    }
  })

  #####################
  # CONTINUOUS VARIABLE FUNCTIONS
  #####################

  value_per_sub_summary <- eventReactive(myReactVals$run_data_button, {
    progress = shiny::Progress$new()
    on.exit(progress$close())
    cont_var_data = get_continuous_var_data()
    if(is.null(cont_var_data)) return(NULL)
    else {
      summary = run_bootstrapping_steps(cont_var_data[converted==T], progress) # only converted users
      value_per_conv_output_table = format_data_summary(summary)
      return(value_per_conv_output_table)
    }
  })

  value_per_event_summary <- eventReactive(myReactVals$run_data_button, {
    progress = shiny::Progress$new()
    on.exit(progress$close())
    cont_var_data = get_continuous_var_data()
    if(is.null(cont_var_data)) return(NULL)
    else {
      summary = run_bootstrapping_steps(cont_var_data[!is.na(user_id)], progress) # every non NA row in the uploaded dataset
      value_per_event_output_table = format_data_summary(summary)
      return(value_per_event_output_table)
    }
  })

  #####################
  # FORMATTING, HEADERS, TEXT
  #####################

  instructions <- reactive({
    # This provides gentle instructions to the user, but note: instructions will not reset if the user iterates new uploads/analysis
    if(input$data_input_type == "upload"){
      if(is.null(input[[testfilepath()]])){
        "Please upload a file."
      } else if(!is.null(input[[testfilepath()]])){
        "Please specify columns"
      }
    }
  })

  instructions_cash <- reactive({
    # This provides gentle instructions to the user, but TODO: instructions will not reset if the user iterates new uploads/analysis
    if(is.null(input[[testfilepath()]])) "Please upload a file."
    else{
      if(is.null(analysis_table())) "Please specify column names, then hit 'Analyze Data'. This may take a few seconds."
      else NULL
    }
  })

  output$summary_text <- renderUI({
    analysis_table = analysis_table()
    if(is.null(analysis_table)){
      NULL
    } else {
      aggs = aggregates_wrapper()
      winning_variant = myReactVals$candidate_winner
      min_effect_size = 1.0 * myReactVals$min_effect_size # as a percent diff from control
      value_remaining = get_value_remaining(aggs, winning_variant, myReactVals$danger_criteria)

      ### SUMMARY
      # calculate prob of win + prob of tie for all variants
      tt = tied_table()
      variant_names = as.character(aggs$arm)
      tied_probability_intermediate= lapply(variant_names, function(x) analysis_table[Variant==x][["Probability of winner (%)"]] + sum(tt[grepl(x, tt[["Tied Winning Variants"]])][[2]]))
      tied_probability = as.data.table(data.frame(Variant = variant_names, prob_win = unlist(tied_probability_intermediate))) # assumes that second colum is the probabilities
      # construct the summary message
      summary_header = h3("Summary:")
      line1 = p('Variant ', strong(winning_variant), 'wins or ties ', strong(paste( round(tied_probability[Variant==winning_variant][["prob_win"]],0),'%', sep = "")), ' of the time given your minimum effect size of ', strong(paste(min_effect_size, "%", sep = "")),".")
      line2 = p(paste("Worst-case scenario (",myReactVals$danger_criteria,"% of the time)", sep=""),
                strong(winning_variant),
                "loses by ",
                strong(paste(value_remaining , "%.", sep = ""))
      )
      summary_message = paste(line1, line2, sep = '<br/>')

      ### RECOMMENDATION
      recommendation_header = h3("Recommendation:")
      recommendation_message = make_recommendation_message(analysis_table, myReactVals$danger_criteria, min_effect_size, tied_probability, aggs)

      return(HTML(paste(summary_header, summary_message, recommendation_header, recommendation_message, sep = '<br/>')))
    }
  })

  headers <- reactive({
    if(is.null(analysis_table()) || ( (input$data_input_type != 'upload') & (input$data_input_type != 'typein') ) ){
      NULL
    } else {
      c("", "Conversion Rate", "Probability distribution", "Differences From Baseline")
    }
  })

  cash_headers <- reactive({
    if(is.null(analysis_table()) || (input$data_input_type != 'upload_cash')){
      NULL
    } else {
      c("", "Subscription Rate", "Probability distribution", "Cash Per Subscription", "Cash Per Event" )
      }
  })

  survey_headers <- reactive({
    if(is.null(value_per_event_summary()) || (input$data_input_type != 'upload_survey')){
      NULL
    } else {
      c("", "Survey Responses Data")
    }
  })
  # define header objects
  output$header2 <- renderText(headers()[2])
  output$header3 <- renderText(headers()[3])
  output$instructions <- renderText(instructions())
  output$header4 <- renderText(headers()[4])
  output$cash_header2 <- renderText(cash_headers()[2])
  output$cash_header3 <- renderText(cash_headers()[3])
  output$cash_header4 <- renderText(cash_headers()[4])
  output$cash_header5 <- renderText(cash_headers()[5])
  output$instructions_cash <- renderText(instructions_cash())
  output$survey_header2 <- renderText(survey_headers()[2])

  output$analysis_table <- renderTable({
    analysis_table()
    }, digits = 2)

  output$analysis_table_cash <- renderTable({
    analysis_table()
  }, digits = 2)

  output$ties_table <- renderTable({
    tied_table()
  }, digits = 2)

  output$group_plot <- renderPlot({
    if(is.null(aggregates_wrapper())) return(invisible())
    else print(make_plot(aggregates_wrapper()))
  })

  output$group_plot_cash <- renderPlot({
    if(is.null(aggregates_wrapper())) return(invisible())
    else print(make_plot(aggregates_wrapper()))
  })

  output$differences_plot <- renderPlot({
    if(is.null(aggregates_wrapper())) return(invisible())
    else print(make_differences_plot(aggregates_wrapper(), myReactVals$baseline_variant))
  })
  output$diff_plot_note <- renderUI({
    if(is.null(analysis_table())){
      NULL
    } else {
      return(p(strong("% Difference")," is defined as  100 * (variant - baseline) / baseline"))
    }
  })

  output$cash_per_sub_table <- renderTable({
    if(myReactVals$doPlot) value_per_sub_summary()
    else NULL
  }, digits = 2)

  output$cash_per_sub_plot <- renderPlot({
    if(is.null(value_per_sub_summary()) || !myReactVals$doPlot) return(invisible())
    else {
      summary_data = as.data.table(value_per_sub_summary())
      summary_data[, up:= `Upper 95% CI`] # this is necessary because make_nice_bar_plot / ggplot doesn't handle column names with spaces well
      summary_data[, lw:= `Lower 95% CI`]
      print(make_nice_bar_plot(summary_data, group_var = "Variant", y = "mean", up = "up", lw = "lw", ylab = "Mean estimate with 95% CI ($)", plot_title = "Cash / Sub Distributions" ))
    }
  })

  output$cash_per_event_table <- renderTable({
    if(myReactVals$doPlot) value_per_event_summary()
    else NULL
  }, digits = 2)

  output$survey_per_event_table <- renderTable({
    if(myReactVals$doPlot) value_per_event_summary()
    else NULL
  }, digits = 2)

  output$cash_per_event_plot <- renderPlot({
    if(is.null(value_per_event_summary())|| !myReactVals$doPlot) return(invisible())
    else {
      summary_data = as.data.table(value_per_event_summary())
      summary_data[, up:= `Upper 95% CI`] # this is necessary because make_nice_bar_plot / ggplot doesn't handle column names with spaces well
      summary_data[, lw:= `Lower 95% CI`]
      print(make_nice_bar_plot(summary_data, group_var = "Variant", y = "mean", up = "up", lw = "lw", ylab = "Mean estimate with 95% CI ($)", plot_title = "Cash / Event Distributions" ))
    }
  })

  output$survey_per_event_plot <- renderPlot({
    if(is.null(value_per_event_summary())|| !myReactVals$doPlot) return(invisible())
    else {
      summary_data = as.data.table(value_per_event_summary())
      summary_data[, up:= `Upper 95% CI`] # this is necessary because make_nice_bar_plot / ggplot doesn't handle column names with spaces well
      summary_data[, lw:= `Lower 95% CI`]
      print(make_nice_bar_plot(summary_data, group_var = "Variant", y = "mean", up = "up", lw = "lw", ylab = "Mean estimate with 95% CI", plot_title = "Survey Response Distributions" ))
    }
  })

}
)
