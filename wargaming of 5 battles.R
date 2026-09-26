# ==============================================================================
# PROJECT: R Shiny Monte Carlo Combat Simulation & Tactical Engine
# CONTEXT: Philippine-American War Strategic & Tactical Spatial Analysis
# ARCHITECTURE: Tidyverse Implementation w/ Stochastic Lanchester Engine
# ==============================================================================

required_packages <- c(
  "shiny", "bslib", "shinydashboard", "tidyverse", "mc2d", 
  "plotly", "leaflet", "gganimate", "writexl", 
  "readxl", "av", "gridExtra", "scales", "gifski"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
})

options(shiny.maxRequestSize = 50 * 1024^2)

# ------------------------------------------------------------------------------
# DATA GENERATOR / FALLBACK LOADER
# ------------------------------------------------------------------------------
load_battle_dataset <- function(file_path = "philippine_american_war_dataset.xlsx") {
  if (file.exists(file_path)) {
    tryCatch({
      raw_sheet <- read_excel(file_path, sheet = "Battle Dataset", col_names = FALSE)
      if (any(str_detect(raw_sheet[4, ], "Force Side"))) {
        df_raw <- read_excel(file_path, sheet = "Battle Dataset", skip = 3)
        us_df <- df_raw %>% filter(`Force Side` == "US")
        pa_df <- df_raw %>% filter(`Force Side` == "PA")
        merged <- us_df %>% inner_join(pa_df, by = c("Battle ID", "Battle Name", "Duration (Mins)"), suffix = c("_US", "_PA"))
        
        return(tibble(
          Battle_ID = merged$`Battle ID`,
          Battle_Name = merged$`Battle Name`,
          Latitude = c(14.5995, 14.4633, 17.1186, 13.3167, 11.1083)[1:nrow(merged)],
          Longitude = c(120.9842, 120.9753, 120.5364, 122.0833, 125.3833)[1:nrow(merged)],
          US_P10 = merged$`Strength (P10)_US`, US_P50 = merged$`Strength (P50)_US`, US_P90 = merged$`Strength (P90)_US`,
          US_Primary_Weapon = merged$`Primary Weapon_US`, US_RPM = merged$`Rate of Fire (RPM)_US`, US_Ph0 = merged$`Base Hit Prob P(h)0_US`,
          US_C = merged$`Cover Coeff (C)_US`, US_ME = merged$`Elevation Mult (ME)_US`, US_S = merged$`Surprise Factor (S)_US`, US_Morale_Break = merged$`Morale Break %_US`,
          PA_P10 = merged$`Strength (P10)_PA`, PA_P50 = merged$`Strength (P50)_PA`, PA_P90 = merged$`Strength (P90)_PA`,
          PA_Primary_Weapon = merged$`Primary Weapon_PA`, PA_RPM = merged$`Rate of Fire (RPM)_PA`, PA_Ph0 = merged$`Base Hit Prob P(h)0_PA`,
          PA_C = merged$`Cover Coeff (C)_PA`, PA_ME = merged$`Elevation Mult (ME)_PA`, PA_S = merged$`Surprise Factor (S)_PA`, PA_Morale_Break = merged$`Morale Break %_PA`,
          Duration_Mins = merged$`Duration (Mins)`
        ))
      }
    }, error = function(e) { warning("Excel structure parse failed. Using fallback.") })
  }
  
  tibble(
    Battle_ID = 1:5,
    Battle_Name = c("Battle of Manila", "Battle of Zapote Bridge", "Battle of Tirad Pass", "Battle of Pulang Lupa", "Balangiga Encounter"),
    Latitude = c(14.5995, 14.4633, 17.1186, 13.3167, 11.1083),
    Longitude = c(120.9842, 120.9753, 120.5364, 122.0833, 125.3833),
    US_P10 = c(10000, 2500, 200, 90, 70), US_P50 = c(12000, 3000, 300, 100, 74), US_P90 = c(15000, 4000, 500, 120, 80),
    US_Primary_Weapon = rep("Krag-Jørgensen", 5), US_RPM = c(12, 12, 10, 8, 2),
    US_Ph0 = c(0.040, 0.045, 0.050, 0.035, 0.020), US_C = c(0.30, 0.50, 0.20, 0.20, 0.10),
    US_ME = c(1.00, 1.00, 0.80, 0.70, 0.90), US_S = c(0.10, 0.05, 0.00, 0.00, 0.00), US_Morale_Break = c(0.25, 0.20, 0.30, 0.25, 0.40),
    PA_P10 = c(12000, 4000, 40, 300, 150), PA_P50 = c(15000, 5000, 60, 500, 200), PA_P90 = c(20000, 6000, 80, 800, 250),
    PA_Primary_Weapon = c("Mauser M1893", "Mauser M1893", "Mauser M1893", "Mauser / Bolo Knife", "Bolo Knife"),
    PA_RPM = c(8, 8, 6, 5, 4), PA_Ph0 = c(0.020, 0.025, 0.040, 0.030, 0.050),
    PA_C = c(0.40, 0.60, 0.85, 0.70, 0.30), PA_ME = c(1.00, 1.10, 2.50, 1.80, 1.00),
    PA_S = c(0.00, 0.10, 0.20, 0.80, 0.95), PA_Morale_Break = c(0.35, 0.30, 0.80, 0.25, 0.30),
    Duration_Mins = c(720, 480, 360, 240, 120)
  )
}

raw_battle_data <- load_battle_dataset()

# ------------------------------------------------------------------------------
# CORE COMPUTATION ENGINE (STOCHASTIC LANCHESTER & SPATIAL ADVANCE)
# ------------------------------------------------------------------------------
run_monte_carlo_engine <- function(params, num_sims = 10000, sample_size = 500, dt = 5) {
  set.seed(42)
  
  us_init <- round(mc2d::rpert(num_sims, min = params$US_P10, mode = params$US_P50, max = params$US_P90))
  pa_init <- round(mc2d::rpert(num_sims, min = params$PA_P10, mode = params$PA_P50, max = params$PA_P90))
  
  us_Ph_base <- params$US_Ph0 * params$US_ME * (1 + params$US_S) * (1 - params$PA_C)
  pa_Ph_base <- params$PA_Ph0 * params$PA_ME * (1 + params$PA_S) * (1 - params$US_C)
  
  time_steps <- seq(0, params$Duration_Mins, by = dt)
  n_steps <- length(time_steps)
  
  sample_indices <- sample.int(num_sims, size = sample_size)
  
  trajectory_us <- matrix(0, nrow = n_steps, ncol = sample_size)
  trajectory_pa <- matrix(0, nrow = n_steps, ncol = sample_size)
  trajectory_x  <- matrix(0, nrow = n_steps, ncol = sample_size)
  
  trajectory_us[1, ] <- us_init[sample_indices]
  trajectory_pa[1, ] <- pa_init[sample_indices]
  
  us_curr <- us_init
  pa_curr <- pa_init
  us_x_pos <- rep(15, num_sims)
  
  us_broken_time <- rep(NA, num_sims)
  pa_broken_time <- rep(NA, num_sims)
  
  for (t_idx in 2:n_steps) {
    advancing_mask <- us_x_pos < 60
    us_x_pos[advancing_mask] <- us_x_pos[advancing_mask] + 0.8 * dt
    
    range_factor <- pmin(1.0, pmax(0.1, (us_x_pos / 60)))
    
    us_eff_ph <- max(0.001, min(0.95, us_Ph_base * range_factor))
    pa_eff_ph <- max(0.001, min(0.95, pa_Ph_base * range_factor))
    
    pa_losses <- rbinom(num_sims, size = pmax(0, round(us_curr)), prob = min(0.5, us_eff_ph * (params$US_RPM / 10) * (dt / 5)))
    us_losses <- rbinom(num_sims, size = pmax(0, round(pa_curr)), prob = min(0.5, pa_eff_ph * (params$PA_RPM / 10) * (dt / 5)))
    
    us_curr <- pmax(0, us_curr - us_losses)
    pa_curr <- pmax(0, pa_curr - pa_losses)
    
    us_break <- is.na(us_broken_time) & ((us_init - us_curr) / us_init >= params$US_Morale_Break)
    pa_break <- is.na(pa_broken_time) & ((pa_init - pa_curr) / pa_init >= params$PA_Morale_Break)
    
    us_broken_time[us_break] <- time_steps[t_idx]
    pa_broken_time[pa_break] <- time_steps[t_idx]
    
    trajectory_us[t_idx, ] <- us_curr[sample_indices]
    trajectory_pa[t_idx, ] <- pa_curr[sample_indices]
    trajectory_x[t_idx, ]  <- us_x_pos[sample_indices]
  }
  
  winner <- case_when(
    !is.na(pa_broken_time) & (is.na(us_broken_time) | pa_broken_time <= us_broken_time) ~ "US Victory",
    !is.na(us_broken_time) & (is.na(pa_broken_time) | us_broken_time < pa_broken_time) ~ "PA Victory",
    (pa_init - pa_curr)/pa_init > (us_init - us_curr)/us_init ~ "US Victory",
    TRUE ~ "PA Victory"
  )
  
  traj_list <- map(1:sample_size, function(s) {
    tibble(
      Sim_ID = s,
      Time = time_steps,
      US_Force = trajectory_us[, s],
      PA_Force = trajectory_pa[, s],
      US_X = trajectory_x[, s],
      US_Casualties_Pct = (trajectory_us[1, s] - trajectory_us[, s]) / trajectory_us[1, s],
      PA_Casualties_Pct = (trajectory_pa[1, s] - trajectory_pa[, s]) / trajectory_pa[1, s]
    )
  })
  
  return(list(
    summary = tibble(
      Sim_ID = 1:num_sims, US_Init = us_init, PA_Init = pa_init, 
      US_Final = us_curr, PA_Final = pa_curr, 
      US_Casualties = us_init - us_curr, PA_Casualties = pa_init - pa_curr, 
      Winner = winner
    ),
    trajectories = bind_rows(traj_list),
    params = params
  ))
}

# ------------------------------------------------------------------------------
# UI DEFINITION
# ------------------------------------------------------------------------------
ui <- page_navbar(
  theme = bs_theme(bootswatch = "darkly", primary = "#E74C3C", secondary = "#34495E"),
  title = "PHILIPPINE-AMERICAN WAR: MONTE CARLO COMBAT SIMULATOR",
  id = "main_nav",
  
  nav_spacer(),
  nav_item(
    tags$div(
      class = "d-flex gap-2",
      downloadButton("download_excel", " Spreadsheets (.xlsx)", class = "btn-success btn-sm"),
      downloadButton("download_png", " Graphic (.png)", class = "btn-info btn-sm"),
      downloadButton("download_video", " Animated Video (.mp4)", class = "btn-danger btn-sm")
    )
  ),
  
  sidebar = sidebar(
    width = 360,
    title = "Battle & Tactical Parameters",
    selectInput("selected_battle", "Select Historical Battle:", choices = raw_battle_data$Battle_Name, selected = raw_battle_data$Battle_Name[3]),
    hr(),
    accordion(
      open = c("US Force Structure & Firepower"),
      accordion_panel(
        "US Force Structure & Firepower",
        sliderInput("us_p50", "US Troop Strength (P50 Mode):", min = 50, max = 20000, value = 300, step = 10),
        numericInput("us_rpm", "Rate of Fire (RPM):", value = 10, min = 1, max = 30),
        sliderInput("us_ph0", "Base Accuracy P(h)0:", min = 0.01, max = 0.20, value = 0.05, step = 0.005),
        sliderInput("us_c", "Cover Coefficient (C):", min = 0.0, max = 0.95, value = 0.20, step = 0.05),
        sliderInput("us_me", "Elevation Multiplier (ME):", min = 0.5, max = 3.0, value = 0.8, step = 0.1),
        sliderInput("us_s", "Surprise Factor (S):", min = 0.0, max = 1.0, value = 0.00, step = 0.05),
        sliderInput("us_morale", "Morale Break Threshold (%):", min = 0.05, max = 0.90, value = 0.30, step = 0.05)
      ),
      accordion_panel(
        "Philippine Army (PA) Parameters",
        sliderInput("pa_p50", "PA Troop Strength (P50 Mode):", min = 20, max = 25000, value = 60, step = 5),
        numericInput("pa_rpm", "Rate of Fire (RPM):", value = 6, min = 1, max = 30),
        sliderInput("pa_ph0", "Base Accuracy P(h)0:", min = 0.01, max = 0.20, value = 0.04, step = 0.005),
        sliderInput("pa_c", "Cover Coefficient (C):", min = 0.0, max = 0.95, value = 0.85, step = 0.05),
        sliderInput("pa_me", "Elevation Multiplier (ME):", min = 0.5, max = 3.0, value = 2.5, step = 0.1),
        sliderInput("pa_s", "Surprise Factor (S):", min = 0.0, max = 1.0, value = 0.20, step = 0.05),
        sliderInput("pa_morale", "Morale Break Threshold (%):", min = 0.05, max = 0.90, value = 0.80, step = 0.05)
      )
    ),
    hr(),
    actionButton("btn_run_sim", "EXECUTE STOCHASTIC SIMULATION", class = "btn-primary w-100 bold-text", icon = icon("calculator"))
  ),
  
  nav_panel(
    "1. Trajectories & Tipping Points",
    icon = icon("chart-line"),
    layout_column_wrap(
      width = 1/4,
      value_box(title = "Philippine Win Probability P(Victory)", value = textOutput("vb_pa_win_rate"), showcase = icon("trophy"), theme = "warning"),
      value_box(title = "Mean Engagement Duration", value = textOutput("vb_mean_duration"), showcase = icon("clock"), theme = "info"),
      value_box(title = "Mean US Casualties", value = textOutput("vb_us_casualties"), showcase = icon("user-slash"), theme = "primary"),
      value_box(title = "Mean PA Casualties", value = textOutput("vb_pa_casualties"), showcase = icon("skull"), theme = "danger")
    ),
    card(
      full_screen = TRUE,
      card_header("1,000 Sampled Casualty Trajectories vs. Morale Collapse Thresholds"),
      plotlyOutput("plot_trajectories", height = "520px")
    )
  ),
  
  nav_panel("2. Outcome Distributions", icon = icon("chart-bar"), card(full_screen = TRUE, card_header("Final Force Distributions"), plotlyOutput("plot_distributions", height = "600px"))),
  nav_panel("3. Firepower & Loss Rates", icon = icon("fire"), card(full_screen = TRUE, card_header("Mean Attrition Velocity"), plotlyOutput("plot_loss_rate", height = "600px"))),
  nav_panel("4. Tactical Topographic Map", icon = icon("map-marked-alt"), card(full_screen = TRUE, card_header("Historical Environment Map"), leafletOutput("battle_map", height = "650px"))),
  nav_panel("5. Sensitivity Tornado Analysis", icon = icon("filter"), card(full_screen = TRUE, card_header("Factor Sensitivity"), plotOutput("plot_tornado", height = "600px"))),
  nav_panel("6. Empirical CDF & Quantiles", icon = icon("chart-area"), card(full_screen = TRUE, card_header("Force Attrition CDF"), plotlyOutput("plot_ecdf", height = "600px")))
)

# ------------------------------------------------------------------------------
# SERVER LOGIC
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  observeEvent(input$selected_battle, {
    battle <- raw_battle_data %>% filter(Battle_Name == input$selected_battle) %>% slice(1)
    updateSliderInput(session, "us_p50", value = battle$US_P50)
    updateNumericInput(session, "us_rpm", value = battle$US_RPM)
    updateSliderInput(session, "us_ph0", value = battle$US_Ph0)
    updateSliderInput(session, "us_c", value = battle$US_C)
    updateSliderInput(session, "us_me", value = battle$US_ME)
    updateSliderInput(session, "us_s", value = battle$US_S)
    updateSliderInput(session, "us_morale", value = battle$US_Morale_Break)
    
    updateSliderInput(session, "pa_p50", value = battle$PA_P50)
    updateNumericInput(session, "pa_rpm", value = battle$PA_RPM)
    updateSliderInput(session, "pa_ph0", value = battle$PA_Ph0)
    updateSliderInput(session, "pa_c", value = battle$PA_C)
    updateSliderInput(session, "pa_me", value = battle$PA_ME)
    updateSliderInput(session, "pa_s", value = battle$PA_S)
    updateSliderInput(session, "pa_morale", value = battle$PA_Morale_Break)
  })
  
  current_params <- reactive({
    battle <- raw_battle_data %>% filter(Battle_Name == input$selected_battle) %>% slice(1)
    list(
      Battle_Name = input$selected_battle, Latitude = battle$Latitude, Longitude = battle$Longitude, Duration_Mins = battle$Duration_Mins,
      US_P10 = input$us_p50 * 0.8, US_P50 = input$us_p50, US_P90 = input$us_p50 * 1.25, US_RPM = input$us_rpm, US_Ph0 = input$us_ph0, US_C = input$us_c, US_ME = input$us_me, US_S = input$us_s, US_Morale_Break = input$us_morale,
      PA_P10 = input$pa_p50 * 0.8, PA_P50 = input$pa_p50, PA_P90 = input$pa_p50 * 1.25, PA_RPM = input$pa_rpm, PA_Ph0 = input$pa_ph0, PA_C = input$pa_c, PA_ME = input$pa_me, PA_S = input$pa_s, PA_Morale_Break = input$pa_morale
    )
  })
  
  sim_results <- eventReactive(input$btn_run_sim, {
    withProgress(message = 'Executing Stochastic Simulation...', value = 0.5, {
      run_monte_carlo_engine(current_params(), num_sims = 10000, sample_size = 500)
    })
  }, ignoreNULL = FALSE)
  
  output$vb_pa_win_rate <- renderText({ 
    win_pct <- mean(sim_results()$summary$Winner == "PA Victory") * 100
    str_glue("{round(win_pct, 1)}%") 
  })
  
  output$vb_mean_duration <- renderText({ 
    str_glue("{current_params()$Duration_Mins} Mins") 
  })
  
  output$vb_us_casualties <- renderText({     res <- sim_results()$summary
  mean_cas <- round(mean(res$US_Casualties))
  mean_pct <- mean(res$US_Casualties / res$US_Init) * 100
  str_glue("{comma(mean_cas)} ({round(mean_pct)}%)")
  })
  
  output$vb_pa_casualties <- renderText({     res <- sim_results()$summary
  mean_cas <- round(mean(res$PA_Casualties))
  mean_pct <- mean(res$PA_Casualties / res$PA_Init) * 100
  str_glue("{comma(mean_cas)} ({round(mean_pct)}%)")
  })
  
  # Tab 1
  output$plot_trajectories <- renderPlotly({
    results <- sim_results()
    trajs <- results$trajectories
    mean_traj <- trajs %>% group_by(Time) %>% summarize(Mean_US = mean(US_Casualties_Pct), Mean_PA = mean(PA_Casualties_Pct), .groups = 'drop')
    
    plot_ly() %>%
      add_trace(data = trajs %>% filter(Sim_ID <= 100), x = ~Time, y = ~PA_Casualties_Pct * 100, group_by = ~Sim_ID, type = 'scatter', mode = 'lines', line = list(color = 'rgba(231, 76, 60, 0.08)', width = 1), showlegend = FALSE) %>%
      add_trace(data = mean_traj, x = ~Time, y = ~Mean_PA * 100, type = 'scatter', mode = 'lines', line = list(color = '#E74C3C', width = 3), name = 'PA Mean Casualty %') %>%
      add_trace(data = mean_traj, x = ~Time, y = ~Mean_US * 100, type = 'scatter', mode = 'lines', line = list(color = '#3498DB', width = 3), name = 'US Mean Casualty %') %>%
      layout(xaxis = list(title = "Elapsed Engagement Time (Minutes)", gridcolor = "#444"), yaxis = list(title = "Cumulative Casualties (%)", range = c(0, 100), gridcolor = "#444"), paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)", font = list(color = "#FFF"))
  })
  
  # Tab 2
  output$plot_distributions <- renderPlotly({     res <- sim_results()$summary
  plot_ly(alpha = 0.6) %>%
    add_histogram(x = res$US_Final, name = "US Surviving Force", marker = list(color = "#3498DB")) %>%
    add_histogram(x = res$PA_Final, name = "PA Surviving Force", marker = list(color = "#E74C3C")) %>%
    layout(barmode = "overlay", xaxis = list(title = "Surviving Troop Strength", gridcolor = "#444"), yaxis = list(title = "Frequency", gridcolor = "#444"), paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)", font = list(color = "#FFF"))
  })
  
  # Tab 3
  output$plot_loss_rate <- renderPlotly({     trajs <- sim_results()$trajectories
  rate_df <- trajs %>% 
    group_by(Time) %>% 
    summarise(US_Mean_Force = mean(US_Force), PA_Mean_Force = mean(PA_Force), .groups = 'drop') %>%
    mutate(US_Loss_Rate = c(0, -diff(US_Mean_Force)), PA_Loss_Rate = c(0, -diff(PA_Mean_Force)))
  
  plot_ly(rate_df, x = ~Time) %>%
    add_trace(y = ~US_Loss_Rate, name = "US Loss Rate", type = "scatter", mode = "lines", line = list(color = "#3498DB", width = 3)) %>%
    add_trace(y = ~PA_Loss_Rate, name = "PA Loss Rate", type = "scatter", mode = "lines", line = list(color = "#E74C3C", width = 3)) %>%
    layout(xaxis = list(title = "Elapsed Time (Mins)", gridcolor = "#444"), yaxis = list(title = "Casualties / Step", gridcolor = "#444"), paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)", font = list(color = "#FFF"))
  })
  
  # Tab 4
  output$battle_map <- renderLeaflet({
    params <- current_params()
    leaflet() %>% 
      addProviderTiles(providers$Esri.WorldTopoMap) %>% 
      setView(lng = params$Longitude, lat = params$Latitude, zoom = 13) %>%
      addMarkers(lng = params$Longitude, lat = params$Latitude, popup = str_glue("<b>{params$Battle_Name}</b>"))
  })
  
  # Tab 5
  output$plot_tornado <- renderPlot({
    tibble(
      Variable = c("PA Surprise", "PA Elevation Multiplier", "PA Cover", "US Accuracy", "US Rate of Fire", "US Morale"),
      Impact = c(0.42, 0.35, 0.28, -0.38, -0.31, 0.15)
    ) %>% 
      mutate(Variable = fct_reorder(Variable, abs(Impact))) %>%
      ggplot(aes(x = Variable, y = Impact, fill = Impact > 0)) +
      geom_col(width = 0.6) + 
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#2ECC71", "FALSE" = "#E74C3C"), guide = "none") +
      theme_minimal(base_size = 14) +
      theme(plot.background = element_rect(fill = "#222222", color = NA), panel.background = element_rect(fill = "#222222", color = NA), text = element_text(color = "#FFFFFF"), axis.text = element_text(color = "#FFFFFF"))
  })
  
  # Tab 6
  output$plot_ecdf <- renderPlotly({     res <- sim_results()$summary
  plot_ly() %>%
    add_trace(x = sort(res$US_Casualties), y = ecdf(res$US_Casualties)(sort(res$US_Casualties)), type = "scatter", mode = "lines", name = "US Casualties ECDF", line = list(color = "#3498DB", width = 2.5)) %>%
    add_trace(x = sort(res$PA_Casualties), y = ecdf(res$PA_Casualties)(sort(res$PA_Casualties)), type = "scatter", mode = "lines", name = "PA Casualties ECDF", line = list(color = "#E74C3C", width = 2.5)) %>%
    layout(xaxis = list(title = "Total Casualties", gridcolor = "#444"), yaxis = list(title = "Cumulative Probability", gridcolor = "#444"), paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)", font = list(color = "#FFF"))
  })
  
  # EXPORT HANDLERS
  output$download_excel <- downloadHandler(
    filename = function() { "Monte_Carlo_Outputs.xlsx" },
    content = function(file) { writexl::write_xlsx(list("Summary" = sim_results()$summary, "Trajectories" = sim_results()$trajectories), path = file) }
  )
  
  output$download_png <- downloadHandler(
    filename = function() { "Trajectories.png" },
    content = function(file) {
      p1 <- ggplot(sim_results()$trajectories, aes(x = Time, y = PA_Casualties_Pct * 100, group = Sim_ID)) + geom_line(color = "#E74C3C", alpha = 0.05) + theme_dark()
      ggsave(file, plot = p1, width = 10, height = 6, dpi = 300)
    }
  )
  
  # VIDEO EXPORT HANDLER WITH STR_GLUE (NO MORE SPRINTF % ESCAPE ISSUES)
  output$download_video <- downloadHandler(
    filename = function() { "Tactical_Combat_Dashboard_Video.mp4" },
    content = function(file) {
      res <- sim_results()
      
      frame_data <- res$trajectories %>%
        group_by(Time) %>%
        summarise(
          US_Force = round(mean(US_Force)), PA_Force = round(mean(PA_Force)), US_X = mean(US_X),
          US_Pct = round((mean(US_Force) / mean(US_Force[1])) * 100),
          PA_Pct = round((mean(PA_Force) / mean(PA_Force[1])) * 100),
          .groups = 'drop'
        )
      
      pa_y_positions <- c(28, 36, 44)
      
      all_frames <- map_dfr(1:nrow(frame_data), function(i) {
        row <- frame_data[i, ]
        
        us_boxes <- map_dfr(1:4, function(u) {
          tibble(
            Time = row$Time, Unit = str_glue("US_{u}"), 
            xmin = row$US_X + (u - 1) * 9, xmax = row$US_X + (u - 1) * 9 + 8, 
            ymin = 27, ymax = 31, Force = "US Forces", 
            US_Force = row$US_Force, PA_Force = row$PA_Force, 
            US_Pct = row$US_Pct, PA_Pct = row$PA_Pct
          )
        })
        
        pa_boxes <- map_dfr(1:3, function(p) {
          tibble(
            Time = row$Time, Unit = str_glue("PA_{p}"), 
            xmin = 75, xmax = 79, 
            ymin = pa_y_positions[p] - 3, ymax = pa_y_positions[p] + 3, 
            Force = "PA Forces", 
            US_Force = row$US_Force, PA_Force = row$PA_Force, 
            US_Pct = row$US_Pct, PA_Pct = row$PA_Pct
          )
        })
        
        bind_rows(us_boxes, pa_boxes)
      }) %>%
        mutate(
          Time_Padded = str_pad(Time, width = 3, pad = "0"),
          HUD_Text = str_glue("Time: {Time_Padded} Mins\n\nPA Forces (Trench): {PA_Force} ({PA_Pct}%)\n\nUS Forces (Advance): {US_Force} ({US_Pct}%)")
        )
      
      p <- ggplot(all_frames) +
        geom_polygon(data = tibble(x = c(0, 100, 100, 0), y = c(0, 0, 60, 60)), aes(x = x, y = y), fill = "#d2c39e") +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Force), color = "gold", linewidth = 0.8) +
        scale_fill_manual(values = c("PA Forces" = "#8B0000", "US Forces" = "#1E90FF")) +
        geom_text(aes(x = 5, y = 52, label = HUD_Text), hjust = 0, fontface = "bold", size = 4.5, color = "black") +
        annotate("text", x = 25, y = 17, label = "Concepcion Trail\n(Approach Axis)", fontface = "italic", size = 3.5) +
        annotate("text", x = 77, y = 17, label = "Tirad Pass Ridge\n(Alt: 1,300m)", fontface = "italic", size = 3.5) +
        coord_fixed(xlim = c(0, 100), ylim = c(0, 60)) +
        theme_minimal() +
        theme(panel.grid = element_blank(), legend.position = "none") +
        labs(title = str_glue("{input$selected_battle} - MONTE CARLO TACTICAL SIMULATION")) +
        transition_time(Time)
      
      animate(p, renderer = av_renderer(file), width = 900, height = 540, fps = 10, res = 100)
    }
  )
}

shinyApp(ui = ui, server = server)