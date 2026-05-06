randomize_blocker_digger_resp <- FALSE

# load data
data_dir <- ""

dvw_files <- list.files(data_dir, pattern = "\\.dvw$", full.names = TRUE)

plays_list <- list()
for (file in dvw_files) {
  print(file)
  data <- try(
    datavolley::dv_read(file, insert_technical_timeouts = FALSE)
  )
  if ("datavolley" %in% class(data)) {
    plays_list[[file]] <- datavolley::plays(data)
  }
}

plays <- dplyr::bind_rows(plays_list) |>
  tibble::as_tibble()

team <- plays |>
  dplyr::filter(!is.na(team_id)) |>
  dplyr::distinct(team_id, team_name = team)

player <- plays |>
  dplyr::filter(!is.na(player_id), player_id != "unknown player") |>
  dplyr::count(player_id, player_name, team_id) |>
  dplyr::arrange(-n) |>
  dplyr::group_by(player_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(player_id, player_name, team_id)

team_lineup_away <- plays |>
  dplyr::select(match_id, set_number, team_id = visiting_team_id, dplyr::starts_with("visiting_player_id")) |>
  dplyr::distinct() |>
  tidyr::pivot_longer(cols = dplyr::starts_with("visiting_player_id"), values_to = "player_id") |>
  dplyr::distinct(match_id, set_number, team_id, player_id) |>
  dplyr::filter(player_id != "unkown player")

team_lineup_home <- plays |>
  dplyr::select(match_id, set_number, team_id = home_team_id, dplyr::starts_with("home_player_id")) |>
  dplyr::distinct() |>
  tidyr::pivot_longer(cols = dplyr::starts_with("home_player_id"), values_to = "player_id") |>
  dplyr::distinct(match_id, set_number, team_id, player_id) |>
  dplyr::filter(player_id != "unkown player")

team_lineup <- dplyr::bind_rows(team_lineup_away, team_lineup_home)

libero <- plays |>
  dplyr::filter(skill %in% c("Serve", "Reception", "Dig")) |>
  dplyr::count(match_id, set_number, team_id, player_id) |>
  dplyr::anti_join(team_lineup, by = c("match_id", "set_number", "team_id", "player_id")) |>
  dplyr::group_by(match_id, set_number, team_id) |>
  dplyr::arrange(-n) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(match_id, set_number, team_id, libero_id = player_id)

# estimate the Markov chain 
#from https://stackoverflow.com/questions/24862046/
cumpaste <- function(x, .sep = " ") {
  concat <- paste(x, collapse = .sep)
  return(substring(concat, 1L, cumsum(c(nchar(x[[1L]]), nchar(x[-1L]) + nchar(.sep)))))
}

contact <- plays |>
  dplyr::mutate(
    team_id_offense = team_id,
    team_id_defense = ifelse(team_id_offense == home_team_id, visiting_team_id, home_team_id),
    serve_receive = ifelse(team == serving_team, "S", "R"),
    is_volley_end = dplyr::coalesce(
      dplyr::lead(point, 1) | (team_touch_id != dplyr::lead(team_touch_id, 1)),
      TRUE
    ),
    abbrev = dplyr::case_when(
      point ~ "P",
      is.na(skill) ~ NA,
      grepl("Unknown", skill) ~ NA,
      skill == "Serve" ~ "SV",
      skill == "Attack" ~ paste0("A", attack_code),
      TRUE ~ paste0(substring(skill, 1, 1), evaluation_code)
    ),
    player_id_lead_1 = dplyr::lead(player_id, 1),
    skill_lead_1 = dplyr::lead(skill, 1),
    evaluation_code_lead_1 = dplyr::lead(evaluation_code, 1),
  ) |>
  dplyr::filter(!is.na(abbrev), !is.na(serving_team)) |>
  dplyr::group_by(match_id, point_id, team_touch_id, point) |>
  # define state of markov chain
  dplyr::mutate(
    state = paste0(serve_receive[1], "_", cumpaste(abbrev, .sep = ".")),
    num_contacts = length(abbrev)  
  ) |>
  dplyr::ungroup()

# data quality checks: can't have more than four touches (including block), and
# serves and attacks must end the volley (may want to relax this assumption for attacks).
# for any points where this is violated, this point is discarded
bad_data <- contact |>
  dplyr::group_by(state) |>
  dplyr::mutate(n = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::filter(
    (num_contacts > 4) |
      (abbrev %in% c("SV", "A") & !is_volley_end) |
      (n < 10)
  ) |>
  dplyr::distinct(match_id, point_id)

# count transitions between each pair of states
state_transition <- contact |>
  dplyr::anti_join(bad_data, by = c("match_id", "point_id")) |>
  dplyr::count(state, next_state = dplyr::lead(state, 1)) |>
  # force points to be terminating states
  dplyr::mutate(
    next_state = dplyr::case_when(
      state == "R_P" ~ "R_P",
      state == "S_P" ~ "S_P",
      TRUE ~ next_state
    )
  )

#convert transition counts into transition probabilities
state_transition_wide <- state_transition |>
  # add a row so that S_SV appears once in the next_state column (for creating transition matrix)
  dplyr::bind_rows(tibble::tibble(state = "S_SV", next_state = "S_SV", n = 0)) |>
  # make sure there are no repeated state-next state pairs
  dplyr::group_by(state, next_state) |>
  dplyr::summarize(n = sum(n), .groups = "drop") |>
  dplyr::group_by(state) |>
  dplyr::transmute(state, next_state, prob = n / sum(n)) |>
  dplyr::arrange(next_state) |>
  tidyr::pivot_wider(names_from = next_state, values_from = prob, values_fill = 0) |>
  dplyr::arrange(state)

# convert dataframe to matrix for multiplication
state_transition_matrix <- state_transition_wide |>
  tibble::column_to_rownames("state") |>
  as.matrix()

state_transition_matrix_limit <- state_transition_matrix

for (i in 1:100) {
  state_transition_matrix_limit <- state_transition_matrix_limit %*% state_transition_matrix
}

#extract terminal state probabilities for each starting state
sideout_prob <- state_transition_matrix_limit[, c("S_P", "R_P")] |>
  tibble::as_tibble() |>
  tibble::add_column(state = rownames(state_transition_matrix), .before = 1) |>
  dplyr::transmute(
    state,
    sideout_prob = R_P
  )

win_prob_from_dig <- sideout_prob |>
  dplyr::filter(nchar(state) == 4, substring(state, 1, 3) == "R_D") |>
  dplyr::transmute(eval = substring(state, 4), win_prob = sideout_prob)

contact_sideout_prob <- contact |>
  dplyr::left_join(sideout_prob, by = "state") |>
  dplyr::mutate(sideout_prob_lead_1 = dplyr::lead(sideout_prob, 1))

# distribute credit/debit to individual players 
serve_data <- contact_sideout_prob |>
  dplyr::filter(
    skill == "Serve",
    skill_lead_1 == "Reception" | evaluation_code == "=",
    !is.na(sideout_prob_lead_1),
    player_id != "unknown player",
    !is.na(home_team_id),
    !is.na(visiting_team_id)
  ) |>
  dplyr::mutate(
    server_id = player_id,
    receiver_id = player_id_lead_1,
    is_error = evaluation_code == "=",
    no_error = !is_error,
    exp_no_error = weighted.mean(sideout_prob_lead_1, w = !is_error),
    pg_is_error = ifelse(is_error, 1, exp_no_error) - sideout_prob,
    pg_no_error = ifelse(is_error, 0, sideout_prob_lead_1 - exp_no_error)
  )

pg_serve <- serve_data |>
  dplyr::select(-player_name, -team_id) |>
  dplyr::left_join(dplyr::select(player, player_id, player_name, team_id), by = "player_id") |>
  dplyr::left_join(team, by = "team_id") |>
  dplyr::group_by(player_id, player_name, team_name) |>
  dplyr::summarize(
    pg = -mean(pg_is_error + pg_no_error),
    n = dplyr::n(),
    pct_is_error = mean(is_error),
    pg_is_error = -mean(pg_is_error),
    n_no_error = sum(no_error),
    pct_ace = weighted.mean(evaluation_code == "#", w = no_error),
    pg_no_error = -mean(pg_no_error),
    .groups = "drop"
  ) |>
  dplyr::arrange(-pg)

pg_reception <- serve_data |>
  dplyr::select(-player_id, -player_name, -team_id) |>
  dplyr::rename(player_id = receiver_id) |>
  dplyr::left_join(dplyr::select(player, player_id, player_name, team_id), by = "player_id") |>
  dplyr::left_join(team, by = "team_id") |>
  dplyr::group_by(player_id, player_name, team_name) |>
  dplyr::summarize(
    n_reception = sum(no_error),
    pct_ace = weighted.mean(evaluation_code == "#", w = no_error),
    pg_reception = weighted.mean(pg_no_error, w = no_error),
    .groups = "drop"
  )

attack <- contact_sideout_prob |>
  dplyr::left_join(libero, by = c("match_id", "set_number", "home_team_id" = "team_id")) |>
  dplyr::rename(home_libero_id = libero_id) |>
  dplyr::left_join(libero, by = c("match_id", "set_number", "visiting_team_id" = "team_id")) |>
  dplyr::rename(visiting_libero_id = libero_id) |>
  dplyr::mutate(
    is_error = evaluation_code == "=",
    is_no_error = !is_error,
    is_block = !is_error & dplyr::coalesce(dplyr::lead(skill, 1) == "Block", FALSE),
    is_no_block = !is_error & !is_block,
    is_block_error = is_block & dplyr::lead(evaluation_code, 1) == "=",
    is_block_no_error = is_block & !is_block_error,
    is_block_return = is_block_no_error & (
      dplyr::lead(evaluation_code, 1) == "#" |
        dplyr::coalesce(dplyr::lead(skill, 2) == "Dig", FALSE) & team_id == dplyr::lead(team_id, 2)
    ),
    is_block_through = is_block_no_error & (
      dplyr::coalesce(dplyr::lead(skill, 2) == "Dig", FALSE) & team_id != dplyr::lead(team_id, 2)
    ),
    outcome_eval = dplyr::case_when(
      # attack error
      is_error ~ "=",
      # clean attack (no block)
      !is_error & !is_block & dplyr::lead(skill, 1) %in% c("Dig", NA) ~ dplyr::coalesce(dplyr::lead(evaluation_code, 1), "="),
      # block error
      is_block_error ~ "=",
      # block touch and through
      is_block & !is_block_error & !is_block_return & dplyr::lead(skill, 2) %in% c("Dig", NA) ~ dplyr::coalesce(dplyr::lead(evaluation_code, 2), "="),
      # block touch and return
      is_block_return & dplyr::lead(skill, 2) %in% c("Dig", NA) ~ dplyr::coalesce(dplyr::lead(evaluation_code, 2), "=")
    ),
    # outcome evaluation code is expressed from which team's perspective?
    outcome_team = ifelse(is_error | is_block_return, "offense", "defense"),
    defense_setter_position = ifelse(team_id == home_team_id, visiting_setter_position, home_setter_position),
    defense_player_id1 = ifelse(team_id == home_team_id, visiting_player_id1, home_player_id1),
    defense_player_id2 = ifelse(team_id == home_team_id, visiting_player_id2, home_player_id2),
    defense_player_id3 = ifelse(team_id == home_team_id, visiting_player_id3, home_player_id3),
    defense_player_id4 = ifelse(team_id == home_team_id, visiting_player_id4, home_player_id4),
    defense_player_id5 = ifelse(team_id == home_team_id, visiting_player_id5, home_player_id5),
    defense_player_id6 = ifelse(team_id == home_team_id, visiting_player_id6, home_player_id6),
    defense_libero_id = ifelse(team_id == home_team_id, visiting_libero_id, home_libero_id),
    offense_setter_position = ifelse(team_id == visiting_team_id, visiting_setter_position, home_setter_position),
    offense_player_id1 = ifelse(team_id == visiting_team_id, visiting_player_id1, home_player_id1),
    offense_player_id2 = ifelse(team_id == visiting_team_id, visiting_player_id2, home_player_id2),
    offense_player_id3 = ifelse(team_id == visiting_team_id, visiting_player_id3, home_player_id3),
    offense_player_id4 = ifelse(team_id == visiting_team_id, visiting_player_id4, home_player_id4),
    offense_player_id5 = ifelse(team_id == visiting_team_id, visiting_player_id5, home_player_id5),
    offense_player_id6 = ifelse(team_id == visiting_team_id, visiting_player_id6, home_player_id6),
    offense_libero_id = ifelse(team_id == visiting_team_id, visiting_libero_id, home_libero_id),
    defense_back_right_player_id = dplyr::case_when(
      defense_setter_position == 1 ~ defense_player_id1,
      defense_setter_position == 2 ~ defense_player_id5,
      defense_setter_position == 3 ~ defense_player_id6,
      defense_setter_position == 4 ~ defense_player_id1,
      defense_setter_position == 5 ~ defense_player_id5,
      defense_setter_position == 6 ~ defense_player_id6
    ),
    defense_front_right_player_id = dplyr::case_when(
      defense_setter_position == 1 ~ defense_player_id4,
      defense_setter_position == 2 ~ defense_player_id2,
      defense_setter_position == 3 ~ defense_player_id3,
      defense_setter_position == 4 ~ defense_player_id4,
      defense_setter_position == 5 ~ defense_player_id2,
      defense_setter_position == 6 ~ defense_player_id3
    ),
    defense_front_middle_player_id = dplyr::case_when(
      defense_setter_position == 1 ~ defense_player_id3,
      defense_setter_position == 2 ~ defense_player_id4,
      defense_setter_position == 3 ~ defense_player_id2,
      defense_setter_position == 4 ~ defense_player_id3,
      defense_setter_position == 5 ~ defense_player_id4,
      defense_setter_position == 6 ~ defense_player_id2
    ),
    defense_front_left_player_id = dplyr::case_when(
      defense_setter_position == 1 ~ defense_player_id2,
      defense_setter_position == 2 ~ defense_player_id3,
      defense_setter_position == 3 ~ defense_player_id4,
      defense_setter_position == 4 ~ defense_player_id2,
      defense_setter_position == 5 ~ defense_player_id3,
      defense_setter_position == 6 ~ defense_player_id4
    ),
    defense_back_left_player_id = defense_libero_id,
    defense_back_middle_player_id = dplyr::case_when(
      defense_setter_position == 1 ~ defense_player_id5,
      defense_setter_position == 2 ~ defense_player_id6,
      defense_setter_position == 3 ~ defense_player_id1,
      defense_setter_position == 4 ~ defense_player_id5,
      defense_setter_position == 5 ~ defense_player_id6,
      defense_setter_position == 6 ~ defense_player_id1
    ),
    offense_back_right_player_id = dplyr::case_when(
      offense_setter_position == 1 ~ offense_player_id1,
      offense_setter_position == 2 ~ offense_player_id5,
      offense_setter_position == 3 ~ offense_player_id6,
      offense_setter_position == 4 ~ offense_player_id1,
      offense_setter_position == 5 ~ offense_player_id5,
      offense_setter_position == 6 ~ offense_player_id6
    ),
    offense_front_right_player_id = dplyr::case_when(
      offense_setter_position == 1 ~ offense_player_id4,
      offense_setter_position == 2 ~ offense_player_id2,
      offense_setter_position == 3 ~ offense_player_id3,
      offense_setter_position == 4 ~ offense_player_id4,
      offense_setter_position == 5 ~ offense_player_id2,
      offense_setter_position == 6 ~ offense_player_id3
    ),
    offense_front_middle_player_id = dplyr::case_when(
      offense_setter_position == 1 ~ offense_player_id3,
      offense_setter_position == 2 ~ offense_player_id4,
      offense_setter_position == 3 ~ offense_player_id2,
      offense_setter_position == 4 ~ offense_player_id3,
      offense_setter_position == 5 ~ offense_player_id4,
      offense_setter_position == 6 ~ offense_player_id2
    ),
    offense_front_left_player_id = dplyr::case_when(
      offense_setter_position == 1 ~ offense_player_id2,
      offense_setter_position == 2 ~ offense_player_id3,
      offense_setter_position == 3 ~ offense_player_id4,
      offense_setter_position == 4 ~ offense_player_id2,
      offense_setter_position == 5 ~ offense_player_id3,
      offense_setter_position == 6 ~ offense_player_id4
    ),
    offense_back_left_player_id = offense_libero_id,
    offense_back_middle_player_id = dplyr::case_when(
      offense_setter_position == 1 ~ offense_player_id5,
      offense_setter_position == 2 ~ offense_player_id6,
      offense_setter_position == 3 ~ offense_player_id1,
      offense_setter_position == 4 ~ offense_player_id5,
      offense_setter_position == 5 ~ offense_player_id6,
      offense_setter_position == 6 ~ offense_player_id1
    ),
    player_rotation = dplyr::case_when(
      player_id == offense_player_id1 ~ 1,
      player_id == offense_player_id2 ~ 2,
      player_id == offense_player_id3 ~ 3,
      player_id == offense_player_id4 ~ 4,
      player_id == offense_player_id5 ~ 5,
      player_id == offense_player_id6 ~ 6
    ),
    player_position = dplyr::case_when(
      player_id == offense_front_left_player_id ~ "OH",
      player_id == offense_front_middle_player_id ~ "MB",
      player_id == offense_front_right_player_id ~ "RS",
      player_id %in% c(offense_back_left_player_id, offense_back_middle_player_id, offense_back_right_player_id) ~ "BR"
    ),
    player_zone = dplyr::case_when(
      player_id == offense_back_left_player_id ~ "BL",
      player_id == offense_back_middle_player_id ~ "BM",
      player_id == offense_back_right_player_id ~ "BR",
      player_id == offense_front_left_player_id ~ "FL",
      player_id == offense_front_middle_player_id ~ "FM",
      player_id == offense_front_right_player_id ~ "FR",
      TRUE ~ "BL" # libero
    ),
    setter_id = ifelse(dplyr::lag(skill, 1) == "Set", dplyr::lag(player_id, 1), NA),
    blocker_id = ifelse(dplyr::lead(skill, 1) == "Block", dplyr::lead(player_id, 1), NA),
    digger_id = dplyr::case_when(
      dplyr::lead(skill, 1) == "Dig" ~ dplyr::lead(player_id, 1),
      dplyr::lead(skill, 1) == "Block" & dplyr::lead(skill, 2) == "Dig" ~ dplyr::lead(player_id, 2),
      TRUE ~ NA
    )
  )

volley_start <- attack |>
  dplyr::group_by(match_id, point_id, team_touch_id) |>
  dplyr::filter(skill != "Block") |>
  dplyr::slice(1) |>
  dplyr::transmute(match_id, point_id, team_touch_id,
                   volley_start = paste0(substring(skill, 1, 1), evaluation_code)
  )

blocker_responsibility <- attack |>
  dplyr::mutate(
    blocker_resp = ifelse(dplyr::lead(skill, 1) == "Block", dplyr::lead(player_zone, 1), NA)
  ) |>
  dplyr::count(attack_code, blocker_resp) |>
  dplyr::filter(!is.na(blocker_resp)) |>
  dplyr::group_by(attack_code) |>
  dplyr::arrange(-n) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(attack_code, blocker_resp)

digger_responsibility <- attack |>
  dplyr::mutate(
    digger_resp = ifelse(dplyr::lead(skill, 1) == "Dig", dplyr::lead(player_zone, 1), NA)
  ) |>
  dplyr::count(attack_code, end_zone, digger_resp) |>
  dplyr::filter(!is.na(digger_resp), !is.na(attack_code)) |>
  dplyr::group_by(attack_code, end_zone) |>
  dplyr::arrange(-n) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(attack_code, end_zone, digger_resp)

attack_data <- attack |>
  dplyr::mutate(
    skill_lead_1 = dplyr::lead(skill, 1),
    skill_lead_2 = dplyr::lead(skill, 2)
  ) |>
  dplyr::filter(
    skill_lead_1 %in% c("Block", "Dig", NA),
    skill_lead_1 %in% c("Dig", NA) | skill_lead_2 %in% c("Dig", NA),
    skill == "Attack",
    !is_block_no_error | xor(is_block_return, is_block_through)
  ) |>
  dplyr::left_join(volley_start, by = c("match_id", "point_id", "team_touch_id")) |>
  dplyr::left_join(blocker_responsibility, by = "attack_code") |>
  dplyr::left_join(digger_responsibility, by = c("attack_code", "end_zone")) |>
  dplyr::left_join(win_prob_from_dig, by = c("outcome_eval" = "eval")) |>
  dplyr::filter(substring(volley_start, 1, 1) %in% c("D", "R", "F")) |>
  dplyr::mutate(
    blocker_resp = ifelse(
      test = rep(randomize_blocker_digger_resp, length = dplyr::n()),
      yes = sample(c("FL", "FM", "FR"), size = dplyr::n(), replace = TRUE),
      no = blocker_resp
    ),
    digger_resp = ifelse(
      test = rep(randomize_blocker_digger_resp, length = dplyr::n()),
      yes = sample(c("BL", "BM", "BR"), size = dplyr::n(), replace = TRUE),
      no = blocker_resp
    ),
    blocker_id = dplyr::case_when(
      !is.na(blocker_id) ~ blocker_id,
      blocker_resp == "FL" ~ defense_front_left_player_id,
      blocker_resp == "FM" ~ defense_front_middle_player_id,
      blocker_resp == "FR" ~ defense_front_right_player_id
    ),
    digger_id = dplyr::case_when(
      !is.na(digger_id) ~ digger_id,
      digger_resp == "FL" ~ defense_front_left_player_id,
      digger_resp == "FM" ~ defense_front_middle_player_id,
      digger_resp == "FR" ~ defense_front_right_player_id,
      digger_resp == "BL" ~ defense_back_left_player_id,
      digger_resp == "BM" ~ defense_back_middle_player_id,
      digger_resp == "BR" ~ defense_back_right_player_id
    ),
    win_prob = ifelse(outcome_team == "offense", win_prob, 1 - win_prob)
  )

exp_win_prob <- attack_data |>
  dplyr::filter(!is.na(win_prob), !is.na(attack_code)) |>
  dplyr::group_by(volley_start, attack_code) |>
  dplyr::summarize(
    n = dplyr::n(),
    exp = mean(win_prob),
    exp_error = 0,
    exp_no_error = weighted.mean(win_prob, w = is_no_error),
    exp_block = weighted.mean(win_prob, w = is_block),
    exp_no_block = weighted.mean(win_prob, w = is_no_block),
    exp_block_error = 1,
    exp_block_no_error = weighted.mean(win_prob, w = is_block_no_error),
    exp_block_return = weighted.mean(win_prob, w = is_block_return),
    exp_block_through = weighted.mean(win_prob, w = is_block_through),
    .groups = "drop"
  ) |>
  dplyr::filter(n >= 10) |>
  dplyr::arrange(-exp)

attack_data_with_pg <- attack_data |>
  dplyr::filter(
    player_id != "unknown player",
    !is.na(home_team_id), !is.na(visiting_team_id), !is.na(win_prob)
  ) |>
  dplyr::inner_join(exp_win_prob, by = c("volley_start", "attack_code")) |>
  dplyr::mutate(
    pg_is_error = dplyr::case_when(is_error ~ exp_error, is_no_error ~ exp_no_error) - exp,
    pg_is_block = dplyr::case_when(is_block ~ exp_block, is_no_block ~ exp_no_block, TRUE ~ exp_no_error) - exp_no_error,
    pg_is_block_error = dplyr::case_when(is_block_error ~ exp_block_error, is_block_no_error ~ exp_block_no_error, TRUE ~ exp_block) - exp_block,
    pg_is_block_return = dplyr::case_when(is_block_return ~ exp_block_return, is_block_through ~ exp_block_through, TRUE ~ exp_block_no_error) - exp_block_no_error,
    pg_no_block = ifelse(is_no_block, win_prob - exp_no_block, 0),
    pg_block_return = ifelse(is_block_return, win_prob - exp_block_return, 0),
    pg_block_through = ifelse(is_block_through, win_prob - exp_block_through, 0)
  )

# Points gained summaries ----
pg_attack <- attack_data_with_pg |>
  dplyr::left_join(team, by = "team_id") |>
  dplyr::group_by(player_id, player_name, team_name, player_position) |>
  dplyr::summarize(
    n = dplyr::n(),
    pct_is_error = mean(is_error),
    pg_is_error = mean(pg_is_error),
    n_no_error = sum(is_no_error),
    pct_is_block = weighted.mean(is_block, w = is_no_error),
    pg_is_block = mean(pg_is_block),
    n_block = sum(is_block),
    pct_is_block_error = weighted.mean(is_block_error, w = is_block),
    pg_is_block_error = mean(pg_is_block_error),
    n_block_no_error = sum(is_block_no_error),
    pct_is_block_through = weighted.mean(is_block_through, w = is_block_no_error),
    pg_is_block_return = mean(pg_is_block_return),
    n_no_block = sum(is_no_block),
    pct_kill_no_block = weighted.mean(win_prob > 0.99, w = is_no_block),
    pg_no_block = mean(pg_no_block),
    n_block_through = sum(is_block_through),
    pct_kill_block_through = weighted.mean(win_prob > 0.99, w = is_block_through),
    pg_block_through = mean(pg_block_through),
    n_block_return = sum(is_block_return),
    pct_stuff = weighted.mean(win_prob < 0.01, w = is_block_return),
    pg_block_return = mean(pg_block_return),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    pg = pg_is_error + pg_is_block + pg_is_block_error + pg_is_block_return + pg_no_block + pg_block_return + pg_block_through
  ) |>
  dplyr::select(player_id, player_name, team_name, player_position, pg, dplyr::everything()) |>
  dplyr::arrange(-pg)

pg_set <- attack_data_with_pg |>
  dplyr::select(-player_id, -player_name) |>
  dplyr::rename(player_id = setter_id) |>
  dplyr::left_join(dplyr::select(player, player_id, player_name), by = "player_id") |>
  dplyr::left_join(team, by = "team_id") |>
  dplyr::group_by(player_id, player_name, team_name) |>
  dplyr::summarize(
    n = dplyr::n(),
    pct_is_error = mean(is_error),
    pg_is_error = mean(pg_is_error),
    n_no_error = sum(is_no_error),
    pct_is_block = weighted.mean(is_block, w = is_no_error),
    pg_is_block = mean(pg_is_block),
    n_block = sum(is_block),
    pct_is_block_error = weighted.mean(is_block_error, w = is_block),
    pg_is_block_error = mean(pg_is_block_error),
    n_block_no_error = sum(is_block_no_error),
    pct_is_block_through = weighted.mean(is_block_through, w = is_block_no_error),
    pg_is_block_return = mean(pg_is_block_return),
    n_no_block = sum(is_no_block),
    pct_kill_no_block = weighted.mean(win_prob > 0.99, w = is_no_block),
    pg_no_block = mean(pg_no_block),
    n_block_through = sum(is_block_through),
    pct_kill_block_through = weighted.mean(win_prob > 0.99, w = is_block_through),
    pg_block_through = mean(pg_block_through),
    n_block_return = sum(is_block_return),
    pct_stuff = weighted.mean(win_prob < 0.01, w = is_block_return),
    pg_block_return = mean(pg_block_return),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    pg = pg_is_error + pg_is_block + pg_is_block_error + pg_is_block_return + pg_no_block + pg_block_return + pg_block_through
  ) |>
  dplyr::select(player_id, player_name, team_name, pg, dplyr::everything()) |>
  dplyr::arrange(-pg)

pg_block <- attack_data_with_pg |>
  dplyr::select(-player_id, -player_name, -team_id) |>
  dplyr::rename(player_id = blocker_id) |>
  dplyr::left_join(player, by = "player_id") |>
  dplyr::left_join(team, by = "team_id") |>
  dplyr::group_by(player_id, player_name, team_name) |>
  dplyr::summarize(
    n = sum(is_no_error),
    pct_is_block = weighted.mean(is_block, w = is_no_error),
    pg_is_block = -weighted.mean(pg_is_block, w = is_no_error),
    n_block = sum(is_block),
    pct_is_block_error = weighted.mean(is_block_error, w = is_block),
    pg_is_block_error = -weighted.mean(pg_is_block_error, w = is_no_error),
    n_block_no_error = sum(is_block_no_error),
    pct_is_block_through = weighted.mean(is_block_through, w = is_block_no_error),
    pg_is_block_return = -weighted.mean(pg_is_block_return, w = is_no_error),
    n_no_block = sum(is_no_block),
    pct_kill_no_block = weighted.mean(win_prob > 0.99, w = is_no_block),
    pg_no_block = -weighted.mean(pg_no_block, w = is_no_error) / 2,
    n_block_through = sum(is_block_through),
    pct_kill_block_through = weighted.mean(win_prob > 0.99, w = is_block_through),
    pg_block_through = -weighted.mean(pg_block_through, w = is_no_error) / 2,
    n_block_return = sum(is_block_return),
    pct_stuff = weighted.mean(win_prob < 0.01, w = is_block_return),
    pg_block_return = -weighted.mean(pg_block_return, w = is_no_error),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    pg = pg_is_block + pg_is_block_error + pg_is_block_return + pg_no_block + pg_block_return + pg_block_through
  ) |>
  dplyr::select(player_id, player_name, team_name, pg, dplyr::everything()) |>
  dplyr::arrange(-pg)

pg_dig <- attack_data_with_pg |>
  dplyr::select(-player_id, -player_name, -team_id) |>
  dplyr::rename(player_id = digger_id) |>
  dplyr::left_join(player, by = "player_id") |>
  dplyr::left_join(team, by = "team_id") |>
  dplyr::group_by(player_id, player_name, team_name) |>
  dplyr::summarize(
    n_dig = sum(is_no_block + is_block_through),
    pct_kill = weighted.mean(win_prob > 0.99, w = is_no_block + is_block_through),
    pg_dig = -weighted.mean(
      x = pg_no_block + pg_block_through,
      w = is_no_block + is_block_through
    ),
    n_no_block = sum(is_no_block),
    pct_kill_no_block = weighted.mean(win_prob > 0.99, w = is_no_block),
    pg_no_block = -weighted.mean(pg_no_block, w = is_no_block),
    n_block_through = sum(is_block_through),
    pct_kill_block_through = weighted.mean(win_prob > 0.99, w = is_block_through),
    pg_block_through = -weighted.mean(pg_block_through, w = is_block_through),
    .groups = "drop"
  ) |>
  dplyr::arrange(-pg_dig)

pg_pass <- pg_reception |>
  dplyr::left_join(pg_dig, by = c("player_id", "player_name", "team_name")) |>
  dplyr::mutate(
    n_reception = dplyr::coalesce(n_reception, 0),
    n_dig = dplyr::coalesce(n_dig, 0),
    n = n_reception + n_dig,
    pg = (n_reception * pg_reception + n_dig * pg_dig) / n
  ) |>
  dplyr::select(player_id, player_name, team_name, n, pg, dplyr::everything()) |>
  dplyr::arrange(-pg)

sets_played <- contact_sideout_prob |>
  dplyr::filter(!is.na(player_id)) |>
  dplyr::count(player_id, match_id, set_number) |>
  dplyr::count(player_id) |>
  dplyr::rename(sets_played = n)

id_columns <- c("player_id", "player_name", "team_name")

pg_overall <- sets_played |>
  dplyr::full_join(
    y = pg_serve |>
      dplyr::mutate(pg = n * pg) |>
      dplyr::select(dplyr::all_of(id_columns), n_serve = n, pg_serve = pg),
    by = "player_id"
  ) |>
  dplyr::full_join(
    y = pg_pass |>
      dplyr::mutate(pg = n * pg) |>
      dplyr::select(dplyr::all_of(id_columns), n_pass = n, pg_pass = pg),
    by = id_columns
  ) |>
  dplyr::full_join(
    y = pg_set |>
      dplyr::mutate(pg = n * pg) |>
      dplyr::select(dplyr::all_of(id_columns), n_set = n, pg_set = pg),
    by = id_columns
  ) |>
  dplyr::full_join(
    y = pg_attack |>
      dplyr::mutate(pg = n * pg) |>
      dplyr::group_by_at(id_columns) |>
      dplyr::summarize(n_attack = sum(n), pg_attack = sum(pg), .groups = "drop"),
    by = id_columns
  ) |>
  dplyr::full_join(
    y = pg_block |>
      dplyr::mutate(pg = n * pg) |>
      dplyr::select(dplyr::all_of(id_columns), n_block = n, pg_block = pg),
    by = id_columns
  ) |>
  dplyr::mutate_at(
    .vars = dplyr::vars(dplyr::matches("^sets_played|^n_|^pg_")),
    .funs = function(x) tidyr::replace_na(x, 0)
  )

pg_overall_per_set <- pg_overall |>
  dplyr::mutate_at(
    .vars = dplyr::vars(dplyr::matches("^n_|^pg_")),
    .funs = function(x, denominator) x / denominator,
    denominator = pg_overall$sets_played
  ) |>
  dplyr::filter(sets_played > 0) |>
  dplyr::mutate(
    n = n_serve + n_pass + n_set + n_attack + n_block,
    pg = pg_serve + pg_pass + pg_set + pg_attack + pg_block,
  ) |>
  dplyr::select(dplyr::all_of(id_columns), sets_played, n, pg, dplyr::everything()) |>
  dplyr::arrange(-pg)

# write outputs 
pg_overall_per_set |>
  dplyr::filter(sets_played >= 15) |>
  write.csv(file = "~/Downloads/pg_overall.csv", row.names = FALSE, na = "")

pg_serve |>
  dplyr::filter(n >= 20) |>
  write.csv(file = "~/Downloads/pg_serve.csv", row.names = FALSE, na = "")

pg_pass |>
  dplyr::filter(n >= 30) |>
  write.csv(file = "~/Downloads/pg_pass.csv", row.names = FALSE, na = "")

pg_set |>
  dplyr::filter(n >= 100) |>
  write.csv(file = "~/Downloads/pg_set.csv", row.names = FALSE, na = "")

pg_attack |>
  dplyr::filter(n >= 30) |>
  write.csv(file = "~/Downloads/pg_attack.csv", row.names = FALSE, na = "")

pg_block |>
  dplyr::filter(n >= 30) |>
  write.csv(file = "~/Downloads/pg_block.csv", row.names = FALSE, na = "")

# visualizations: Points Gained per Opportunity (all players) 
harvard_blue <- rgb(0.004, 0.149, 0.314, alpha = 0.6, maxColorValue = 1)
harvard_blue_border <- rgb(0.004, 0.149, 0.314, maxColorValue = 1)

pg_hist <- function(x, n, main = "") {
  # expand to one row per contact
  expanded <- rep(x, times = round(n))
  hist(
    x = expanded,
    breaks = c(-1, seq(from = -0.1, to = 0.1, by = 0.01), 1),
    xlim = c(-0.1, 0.1),
    xlab = "",
    ylab = "",
    axes = FALSE,
    col = harvard_blue,
    border = harvard_blue_border,
    main = main
  )
  axis(
    side = 1,
    at = seq(from = -0.1, to = 0.1, by = 0.05),
    labels = c("-10%", "-5%", "0%", "5%", "10%")
  )
  axis(side = 2)
}

{
  pdf("~/Downloads/points_gained_per_opportunity_all.pdf", height = 6, width = 9)
  par(mfrow = c(2, 3))
  par(mar = c(4.1, 4.1, 3.1, 1.1))
  
  pg_serve |>
    dplyr::filter(!is.nan(pg), !is.na(pg)) |>
    with(pg_hist(pg, n, main = "Serve"))
  
  pg_set |>
    dplyr::filter(!is.nan(pg), !is.na(pg)) |>
    with(pg_hist(pg, n, main = "Set"))
  
  pg_attack |>
    dplyr::group_by(player_id, player_name, team_name) |>
    dplyr::summarize(pg = mean(pg, na.rm = TRUE), n = sum(n), .groups = "drop") |>
    dplyr::filter(!is.nan(pg), !is.na(pg)) |>
    with(pg_hist(pg, n, main = "Attack"))
  
  pg_pass |>
    dplyr::filter(!is.nan(pg_reception), !is.na(pg_reception)) |>
    with(pg_hist(pg_reception, n_reception, main = "Reception"))
  
  pg_pass |>
    dplyr::filter(!is.nan(pg_dig), !is.na(pg_dig)) |>
    with(pg_hist(pg_dig, n_dig, main = "Dig"))
  
  pg_block |>
    dplyr::filter(!is.nan(pg), !is.na(pg)) |>
    with(pg_hist(pg, n, main = "Block"))
  
  title(
    xlab = "Points Gained per contact",
    ylab = "Number of contacts",
    outer = TRUE, line = -1
  )
  
  dev.off()
}

#bar charts with player labels (all players, Harvard highlighted) 
pg_bar <- function(data, pg_col, name_col = "player_name", main = "", min_n = NULL, n_col = "n") {
  data <- data |>
    dplyr::filter(!is.nan(!!rlang::sym(pg_col)), !is.na(!!rlang::sym(pg_col))) |>
    dplyr::arrange(!!rlang::sym(pg_col))
  
  if (!is.null(min_n)) {
    data <- data |> dplyr::filter(!!rlang::sym(n_col) >= min_n)
  }
  
  values  <- data[[pg_col]]
  names   <- data[[name_col]]
  is_harv <- data[["team_name"]] == "Harvard University"
  colors  <- dplyr::case_when(
    is_harv & values >= 0  ~ harvard_blue_border,
    is_harv & values < 0   ~ rgb(0.7, 0.1, 0.1),
    !is_harv & values >= 0 ~ rgb(0.6, 0.7, 0.8),
    TRUE                   ~ rgb(0.8, 0.6, 0.6)
  )
  
  par(mar = c(5.1, 10.1, 4.1, 2.1))
  barplot(
    values,
    names.arg = names,
    horiz = TRUE,
    las = 1,
    col = colors,
    border = colors,
    main = main,
    xlab = "Points Gained per Opportunity",
    xlim = c(
      min(values, -0.05) * 1.2,
      max(values, 0.05) * 1.2
    ),
    axes = FALSE,
    cex.names = 0.7
  )
  abline(v = 0, lty = 2, col = "gray40")
  axis(
    side = 1,
    at = seq(-0.15, 0.15, by = 0.05),
    labels = paste0(seq(-15, 15, by = 5), "%")
  )
  legend("bottomright", bty = "n", cex = 0.7,
         legend = c("Harvard", "Opponent"),
         fill = c(harvard_blue_border, rgb(0.6, 0.7, 0.8)),
         border = c(harvard_blue_border, rgb(0.6, 0.7, 0.8)))
}

{
  pdf("~/Downloads/player_pg_barcharts_all.pdf", height = 10, width = 9)
  
  pg_bar(pg_serve, "pg", min_n = 20, main = "Serve: Points Gained per Opportunity")
  pg_bar(pg_set, "pg", min_n = 30, main = "Set: Points Gained per Opportunity")
  pg_bar(pg_attack |> dplyr::group_by(player_id, player_name, team_name) |>
           dplyr::summarize(pg = mean(pg), n = sum(n), .groups = "drop"),
         "pg", min_n = 30, main = "Attack: Points Gained per Opportunity")
  pg_bar(pg_pass |> dplyr::mutate(n = n_reception), "pg_reception",
         min_n = 20, n_col = "n_reception", main = "Reception: Points Gained per Opportunity")
  pg_bar(pg_pass |> dplyr::mutate(n = n_dig), "pg_dig",
         min_n = 20, n_col = "n_dig", main = "Dig: Points Gained per Opportunity")
  pg_bar(pg_block, "pg", min_n = 30, main = "Block: Points Gained per Opportunity")
  pg_bar(pg_overall_per_set |> dplyr::filter(sets_played >= 15), "pg",
         main = "Overall: Points Gained per Set")
  
  dev.off()
}

# mixed effects models for variance decomposition 
attack_data_with_pg <- attack_data_with_pg |>
  dplyr::group_by(player_id) |>
  dplyr::mutate(count = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::mutate(player_id_model = ifelse(count < 50, "0", player_id)) |>
  dplyr::group_by(setter_id) |>
  dplyr::mutate(count = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::mutate(setter_id_model = ifelse(count < 100, "0", setter_id)) |>
  dplyr::group_by(blocker_id) |>
  dplyr::mutate(count = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::mutate(blocker_id_model = ifelse(count < 50, "0", blocker_id)) |>
  dplyr::group_by(digger_id) |>
  dplyr::mutate(count = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::mutate(digger_id_model = ifelse(count < 30, "0", digger_id))

fit_attack_models <- function(data) {
  models <- list()
  models$is_error <- lme4::lmer(
    pg_is_error ~ (1 | player_id_model) + (1 | setter_id_model),
    data = data,
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
  models$is_block <- lme4::lmer(
    pg_is_block ~ (1 | player_id_model) + (1 | setter_id_model) + (1 | blocker_id_model),
    data = data |> dplyr::filter(is_no_error),
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
  models$is_block_error <- lme4::lmer(
    pg_is_block_error ~ (1 | player_id_model) + (1 | setter_id_model) + (1 | blocker_id_model),
    data = data |> dplyr::filter(is_block),
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
  models$is_block_return <- lme4::lmer(
    pg_is_block_return ~ (1 | player_id_model) + (1 | setter_id_model) + (1 | blocker_id_model),
    data = data |> dplyr::filter(is_block_no_error),
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
  models$no_block <- lme4::lmer(
    pg_no_block ~ (1 | player_id_model) + (1 | setter_id_model) + (1 | blocker_id_model) + (1 | digger_id_model),
    data = data |> dplyr::filter(is_no_block),
    control = lme4::lmerControl(calc.derivs = FALSE, optimizer = "Nelder_Mead")
  )
  models$block_through <- lme4::lmer(
    pg_block_through ~ (1 | player_id_model) + (1 | setter_id_model) + (1 | blocker_id_model) + (1 | digger_id_model),
    data = data |> dplyr::filter(is_block_through),
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
  models$block_return <- lme4::lmer(
    pg_block_return ~ (1 | player_id_model) + (1 | setter_id_model) + (1 | blocker_id_model),
    data = data |> dplyr::filter(is_block_return),
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
  return(models)
}

extract_variance_split <- function(models) {
  get_sd <- function(model, term) {
    vc <- lme4::VarCorr(model)
    if (term %in% names(vc)) attr(vc[[term]], "stddev") else 0
  }
  tibble::tibble(
    model        = names(models),
    sd_attacker  = sapply(models, get_sd, term = "player_id_model"),
    sd_setter    = sapply(models, get_sd, term = "setter_id_model"),
    sd_blocker   = sapply(models, get_sd, term = "blocker_id_model"),
    sd_digger    = sapply(models, get_sd, term = "digger_id_model")
  ) |>
    dplyr::mutate(
      total_as     = sd_attacker + sd_setter,
      pct_attacker = sd_attacker / total_as,
      pct_setter   = sd_setter   / total_as,
      total_bd     = sd_blocker + sd_digger,
      pct_blocker  = ifelse(total_bd > 0, sd_blocker / total_bd, NA),
      pct_digger   = ifelse(total_bd > 0, sd_digger  / total_bd, NA)
    )
}

cat("Fitting attack models...\n")
attack_models  <- fit_attack_models(attack_data_with_pg)
variance_split <- extract_variance_split(attack_models)
cat("Done fitting models.\n")

# bootstrap standard errors 
cat("Running 200 bootstrap resamples...\n")
set.seed(42)
match_ids <- unique(plays$match_id)
n_boot    <- 200
boot_results <- vector("list", n_boot)

for (b in 1:n_boot) {
  if (b %% 20 == 0) cat(sprintf("  Bootstrap %d / %d\n", b, n_boot))
  
  sampled_matches <- sample(match_ids, length(match_ids), replace = TRUE)
  
  boot_plays <- purrr::map_dfr(seq_along(sampled_matches), function(i) {
    plays |>
      dplyr::filter(match_id == sampled_matches[i]) |>
      dplyr::mutate(match_id = paste0(match_id, "_", i))
  })
  
  boot_contact <- tryCatch({
    boot_plays |>
      dplyr::mutate(
        serve_receive = ifelse(team == serving_team, "S", "R"),
        is_volley_end = dplyr::coalesce(
          dplyr::lead(point, 1) | (team_touch_id != dplyr::lead(team_touch_id, 1)), TRUE),
        abbrev = dplyr::case_when(
          point ~ "P",
          is.na(skill) ~ NA,
          grepl("Unknown", skill) ~ NA,
          skill == "Serve" ~ "SV",
          skill == "Attack" ~ paste0("A", attack_code),
          TRUE ~ paste0(substring(skill, 1, 1), evaluation_code)
        ),
        player_id_lead_1 = dplyr::lead(player_id, 1),
        skill_lead_1 = dplyr::lead(skill, 1)
      ) |>
      dplyr::filter(!is.na(abbrev), !is.na(serving_team)) |>
      dplyr::group_by(match_id, point_id, team_touch_id, point) |>
      dplyr::mutate(
        state = paste0(serve_receive[1], "_", cumpaste(abbrev, .sep = ".")),
        num_contacts = length(abbrev)
      ) |>
      dplyr::ungroup()
  }, error = function(e) NULL)
  
  if (is.null(boot_contact)) next
  
  boot_bad <- boot_contact |>
    dplyr::group_by(state) |> dplyr::mutate(n = dplyr::n()) |> dplyr::ungroup() |>
    dplyr::filter((num_contacts > 4) | (abbrev %in% c("SV", "A") & !is_volley_end) | (n < 10)) |>
    dplyr::distinct(match_id, point_id)
  
  boot_stw <- tryCatch({
    boot_contact |>
      dplyr::anti_join(boot_bad, by = c("match_id", "point_id")) |>
      dplyr::count(state, next_state = dplyr::lead(state, 1)) |>
      dplyr::mutate(next_state = dplyr::case_when(
        state == "R_P" ~ "R_P", state == "S_P" ~ "S_P", TRUE ~ next_state)) |>
      dplyr::bind_rows(tibble::tibble(state = "S_SV", next_state = "S_SV", n = 0)) |>
      dplyr::group_by(state, next_state) |>
      dplyr::summarize(n = sum(n), .groups = "drop") |>
      dplyr::group_by(state) |>
      dplyr::transmute(state, next_state, prob = n / sum(n)) |>
      dplyr::arrange(next_state) |>
      tidyr::pivot_wider(names_from = next_state, values_from = prob, values_fill = 0) |>
      dplyr::arrange(state)
  }, error = function(e) NULL)
  
  if (is.null(boot_stw)) next
  
  boot_mat <- boot_stw |> tibble::column_to_rownames("state") |> as.matrix()
  boot_mat_lim <- boot_mat
  for (i in 1:100) boot_mat_lim <- boot_mat_lim %*% boot_mat
  
  if (!all(c("S_P", "R_P") %in% colnames(boot_mat_lim))) next
  
  boot_sp <- boot_mat_lim[, c("S_P", "R_P")] |>
    tibble::as_tibble() |>
    tibble::add_column(state = rownames(boot_mat_lim), .before = 1) |>
    dplyr::transmute(state, sideout_prob = R_P)
  
  boot_wfd <- boot_sp |>
    dplyr::filter(nchar(state) == 4, substring(state, 1, 3) == "R_D") |>
    dplyr::transmute(eval = substring(state, 4), win_prob = sideout_prob)
  
  boot_pg <- tryCatch({
    boot_contact |>
      dplyr::left_join(boot_sp, by = "state") |>
      dplyr::mutate(sideout_prob_lead_1 = dplyr::lead(sideout_prob, 1)) |>
      dplyr::filter(skill == "Attack", !is.na(sideout_prob), player_id != "unknown player") |>
      dplyr::mutate(
        is_error        = evaluation_code == "=",
        is_no_error     = !is_error,
        is_block        = !is_error & dplyr::coalesce(dplyr::lead(skill_lead_1 == "Block"), FALSE),
        is_no_block     = !is_error & !is_block,
        is_block_error  = is_block & dplyr::lead(evaluation_code, 1) == "=",
        is_block_no_error = is_block & !is_block_error,
        is_block_return = is_block_no_error & (dplyr::lead(evaluation_code, 1) == "#" |
                                                 dplyr::coalesce(dplyr::lead(skill_lead_1 == "Dig"), FALSE)),
        is_block_through = is_block_no_error & !is_block_return,
        outcome_team    = ifelse(is_error | is_block_return, "offense", "defense"),
        win_prob        = ifelse(outcome_team == "offense", sideout_prob, 1 - sideout_prob)
      ) |>
      dplyr::filter(!is.na(attack_code)) |>
      dplyr::group_by(attack_code) |>
      dplyr::mutate(
        n_ac        = dplyr::n(),
        exp         = mean(win_prob),
        exp_error   = 0,
        exp_no_error = weighted.mean(win_prob, w = is_no_error),
        exp_block   = weighted.mean(win_prob, w = is_block),
        exp_no_block = weighted.mean(win_prob, w = is_no_block),
        exp_block_error = 1,
        exp_block_no_error = weighted.mean(win_prob, w = is_block_no_error),
        exp_block_return = weighted.mean(win_prob, w = is_block_return),
        exp_block_through = weighted.mean(win_prob, w = is_block_through)
      ) |>
      dplyr::ungroup() |>
      dplyr::filter(n_ac >= 10) |>
      dplyr::mutate(
        pg_is_error      = dplyr::case_when(is_error ~ exp_error, is_no_error ~ exp_no_error) - exp,
        pg_is_block      = dplyr::case_when(is_block ~ exp_block, is_no_block ~ exp_no_block, TRUE ~ exp_no_error) - exp_no_error,
        pg_is_block_error = dplyr::case_when(is_block_error ~ exp_block_error, is_block_no_error ~ exp_block_no_error, TRUE ~ exp_block) - exp_block,
        pg_is_block_return = dplyr::case_when(is_block_return ~ exp_block_return, is_block_through ~ exp_block_through, TRUE ~ exp_block_no_error) - exp_block_no_error,
        pg_no_block      = ifelse(is_no_block, win_prob - exp_no_block, 0),
        pg_block_return  = ifelse(is_block_return, win_prob - exp_block_return, 0),
        pg_block_through = ifelse(is_block_through, win_prob - exp_block_through, 0)
      ) |>
      dplyr::group_by(player_id) |>
      dplyr::mutate(count = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        player_id_model  = ifelse(count < 50, "0", player_id),
        setter_id_model  = ifelse(count < 50, "0", player_id),
        blocker_id_model = ifelse(count < 50, "0", player_id),
        digger_id_model  = ifelse(count < 30, "0", player_id)
      )
  }, error = function(e) NULL)
  
  if (is.null(boot_pg) || nrow(boot_pg) == 0) next
  
  boot_models <- tryCatch(suppressWarnings(fit_attack_models(boot_pg)), error = function(e) NULL)
  if (is.null(boot_models)) next
  
  boot_results[[b]] <- tryCatch(extract_variance_split(boot_models), error = function(e) NULL)
}

cat("Bootstrap complete.\n")

boot_df <- dplyr::bind_rows(purrr::compact(boot_results))

boot_se <- boot_df |>
  dplyr::group_by(model) |>
  dplyr::summarize(
    se_attacker = sd(pct_attacker, na.rm = TRUE),
    se_setter   = sd(pct_setter,   na.rm = TRUE),
    se_blocker  = sd(pct_blocker,  na.rm = TRUE),
    se_digger   = sd(pct_digger,   na.rm = TRUE),
    .groups = "drop"
  )

# build and save the table as PNG 
model_labels <- c(
  "is_error"        = "(1) Attack error",
  "is_block"        = "(2) Block indicator",
  "is_block_error"  = "(3) Block error",
  "is_block_return" = "(4) Block-return indicator",
  "no_block"        = "(5) Clean attack outcome",
  "block_through"   = "(6) Block-through outcome",
  "block_return"    = "(7) Block-return outcome"
)

fmt_pct <- function(x) ifelse(is.na(x), "\u2013", paste0(round(x * 100), "%"))
fmt_se  <- function(x) ifelse(is.na(x), "\u2013", paste0("\u00b1", round(x * 100, 1), "%"))

table_data <- variance_split |>
  dplyr::left_join(boot_se, by = "model") |>
  dplyr::mutate(model_label = model_labels[model]) |>
  dplyr::select(model_label, pct_attacker, pct_setter, se_attacker,
                pct_blocker, pct_digger, se_blocker) |>
  dplyr::mutate(
    pct_attacker = fmt_pct(pct_attacker),
    pct_setter   = fmt_pct(pct_setter),
    se_attacker  = fmt_se(se_attacker),
    pct_blocker  = fmt_pct(pct_blocker),
    pct_digger   = fmt_pct(pct_digger),
    se_blocker   = fmt_se(se_blocker)
  )

# Transpose: models become columns, roles become rows
model_short <- c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)")

table_wide <- tibble::tibble(
  Role = c("Attacker", "Setter", "Std err", "Blocker", "Digger", "Std err")
)

for (i in seq_len(nrow(table_data))) {
  col_vals <- c(
    table_data$pct_attacker[i],
    table_data$pct_setter[i],
    table_data$se_attacker[i],
    table_data$pct_blocker[i],
    table_data$pct_digger[i],
    table_data$se_blocker[i]
  )
  table_wide[[model_short[i]]] <- col_vals
}

gt_table <- table_wide |>
  gt::gt() |>
  gt::cols_label(Role = "") |>
  gt::tab_spanner(
    label   = "Attacker vs Setter",
    columns = c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)")
  ) |>
  gt::tab_row_group(label = "Blocker vs Digger", rows = 4:6) |>
  gt::tab_row_group(label = "Attacker vs Setter", rows = 1:3) |>
  gt::row_group_order(groups = c("Attacker vs Setter", "Blocker vs Digger")) |>
  gt::tab_header(
    title    = "Division of Points Gained Between Teammates",
    subtitle = "Harvard volleyball \u2014 attack outcome tree splits"
  ) |>
  gt::tab_style(
    style = gt::cell_text(color = "gray50", style = "italic"),
    locations = gt::cells_body(rows = Role == "Std err")
  ) |>
  gt::opt_table_font(font = gt::google_font("Source Sans Pro")) |>
  gt::tab_options(
    table.background.color    = "#FFFFF8",
    column_labels.font.weight = "bold",
    row_group.font.weight     = "bold",
    heading.title.font.size   = 14
  )

gt::gtsave(gt_table, filename = "~/Downloads/variance_decomposition_harvard.png", zoom = 2)
cat("Table saved to ~/Downloads/variance_decomposition_harvard.png\n")

# scatter plots: rate stat vs PG (all players, Harvard labeled + top/bottom opponents) 
gray_pt  <- rgb(0.6, 0.6, 0.6, alpha = 0.5)
harv_col <- harvard_blue_border

pg_scatter <- function(all_data, harv_data,
                       x_col, y_col,
                       xlab, ylab, main,
                       x_pct = TRUE, y_pct = TRUE,
                       min_n_col = "n", min_n = 20,
                       n_label_opponents = 3) {
  
  all_data <- all_data |>
    dplyr::filter(!!rlang::sym(min_n_col) >= min_n) |>
    dplyr::filter(!is.nan(!!rlang::sym(y_col)), !is.na(!!rlang::sym(y_col)),
                  !is.nan(!!rlang::sym(x_col)), !is.na(!!rlang::sym(x_col)))
  
  harv_data <- all_data |>
    dplyr::filter(team_name == "Harvard University")
  
  opp_data <- all_data |>
    dplyr::filter(team_name != "Harvard University")
  
  # top and bottom opponents by y_col
  opp_label <- dplyr::bind_rows(
    opp_data |> dplyr::arrange(dplyr::desc(!!rlang::sym(y_col))) |> head(n_label_opponents),
    opp_data |> dplyr::arrange(!!rlang::sym(y_col)) |> head(n_label_opponents)
  ) |> dplyr::distinct()
  
  x_all  <- all_data[[x_col]]
  y_all  <- all_data[[y_col]]
  x_harv <- harv_data[[x_col]]
  y_harv <- harv_data[[y_col]]
  x_opp_label <- opp_label[[x_col]]
  y_opp_label <- opp_label[[y_col]]
  
  x_range <- range(x_all, na.rm = TRUE)
  y_range <- range(y_all, na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.08
  y_pad   <- diff(y_range) * 0.15
  
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  plot(
    x_all, y_all,
    col  = gray_pt, pch = 1, cex = 0.8,
    xlim = x_range + c(-x_pad, x_pad),
    ylim = y_range + c(-y_pad, y_pad),
    xlab = xlab, ylab = ylab, main = main,
    axes = FALSE, type = "n"
  )
  
  abline(h = 0, lty = 2, col = "gray60")
  abline(v = mean(x_all, na.rm = TRUE), lty = 2, col = "gray60")
  
  # all opponents as gray circles
  points(opp_data[[x_col]], opp_data[[y_col]], col = gray_pt, pch = 1, cex = 0.8)
  
  # harvard as navy diamonds
  points(x_harv, y_harv, col = harv_col, pch = 18, cex = 1.8)
  
  # label all Harvard players
  for (i in seq_along(x_harv)) {
    text(x_harv[i], y_harv[i], labels = harv_data[["player_name"]][i],
         pos = 4, cex = 0.65, col = harv_col, offset = 0.4)
  }
  
  # label top/bottom opponents in gray
  for (i in seq_along(x_opp_label)) {
    text(x_opp_label[i], y_opp_label[i], labels = opp_label[["player_name"]][i],
         pos = 4, cex = 0.55, col = "gray40", offset = 0.4)
  }
  
  x_at <- pretty(x_all, n = 5)
  y_at <- pretty(y_all, n = 5)
  axis(1, at = x_at,
       labels = if (x_pct) paste0(round(x_at * 100), "%") else round(x_at, 3))
  axis(2, at = y_at,
       labels = if (y_pct) paste0(round(y_at * 100), "%") else round(y_at, 3))
  
  legend("topright", bty = "n", cex = 0.7,
         legend = c("Harvard", "Opponent"),
         pch = c(18, 1),
         col = c(harv_col, gray_pt))
}

{
  pdf("~/Downloads/pg_scatter_harvard.pdf", height = 6, width = 7)
  
  # serve: error rate vs PG
  pg_scatter(
    all_data  = pg_serve,
    harv_data = pg_serve,
    x_col = "pct_is_error", y_col = "pg",
    xlab  = "Serve Error Rate (min. 20 serves)",
    ylab  = "Points Gained per Opportunity",
    main  = "Serve",
    min_n_col = "n", min_n = 20
  )
  
  # reception: ace rate faced vs PG
  pg_scatter(
    all_data  = pg_reception,
    harv_data = pg_reception,
    x_col = "pct_ace", y_col = "pg_reception",
    xlab  = "Ace Rate Faced (min. 20 receptions)",
    ylab  = "Points Gained per Opportunity",
    main  = "Reception",
    min_n_col = "n_reception", min_n = 20
  )
  
  # attack: kill rate vs PG
  pg_attack_combined <- pg_attack |>
    dplyr::group_by(player_id, player_name, team_name) |>
    dplyr::summarize(
      n          = sum(n),
      n_no_block = sum(n_no_block),
      n_kill     = sum(pct_kill_no_block * n_no_block),
      pct_kill   = n_kill / n_no_block,
      pg         = mean(pg),
      .groups = "drop"
    )
  
  pg_scatter(
    all_data  = pg_attack_combined,
    harv_data = pg_attack_combined,
    x_col = "pct_kill", y_col = "pg",
    xlab  = "Kill Rate - No Block (min. 30 attacks)",
    ylab  = "Points Gained per Opportunity",
    main  = "Attack: Efficiency",
    min_n_col = "n", min_n = 30
  )
  
  # attack: opportunities vs PG
  pg_scatter(
    all_data  = pg_attack_combined,
    harv_data = pg_attack_combined,
    x_col = "n", y_col = "pg",
    xlab  = "Attack Opportunities (min. 30 attacks)",
    ylab  = "Points Gained per Opportunity",
    main  = "Attack: Volume vs Quality",
    x_pct = FALSE,
    min_n_col = "n", min_n = 30
  )
  
  # dig: digs per opportunity vs PG 
  pg_pass_dig <- pg_pass |>
    dplyr::filter(!is.na(pg_dig), n_dig >= 20) |>
    dplyr::mutate(digs_per_opp = 1 - pct_kill)
  
  highlight_players <- c("Zoe Leonard", "Paris Winkler")
  
  pg_pass_dig_all       <- pg_pass_dig
  pg_pass_dig_highlight <- pg_pass_dig |>
    dplyr::filter(player_name %in% highlight_players)
  pg_pass_dig_other     <- pg_pass_dig |>
    dplyr::filter(!player_name %in% highlight_players)
  
  {
    x_all   <- pg_pass_dig_all$digs_per_opp
    y_all   <- pg_pass_dig_all$pg_dig
    x_hi    <- pg_pass_dig_highlight$digs_per_opp
    y_hi    <- pg_pass_dig_highlight$pg_dig
    names   <- pg_pass_dig_highlight$player_name
    
    x_range <- range(x_all, na.rm = TRUE)
    y_range <- range(y_all, na.rm = TRUE)
    x_pad   <- diff(x_range) * 0.08
    y_pad   <- diff(y_range) * 0.15
    
    par(mar = c(5.1, 5.1, 4.1, 2.1))
    plot(
      x_all, y_all,
      col  = gray_pt, pch = 1, cex = 0.8,
      xlim = x_range + c(-x_pad, x_pad),
      ylim = y_range + c(-y_pad, y_pad),
      xlab = "Digs per Opportunity (min. 20 opportunities)",
      ylab = "Points Gained per Opportunity",
      main = "Dig",
      axes = FALSE, type = "n"
    )
    abline(h = 0, lty = 2, col = "gray60")
    abline(v = mean(x_all, na.rm = TRUE), lty = 2, col = "gray60")
    points(pg_pass_dig_other$digs_per_opp, pg_pass_dig_other$pg_dig,
           col = gray_pt, pch = 1, cex = 0.8)
    points(x_hi, y_hi, col = harv_col, pch = 18, cex = 2)
    for (i in seq_along(x_hi)) {
      text(x_hi[i], y_hi[i], labels = names[i],
           pos = 4, cex = 0.75, col = harv_col, offset = 0.4)
    }
    x_at <- pretty(x_all, n = 5)
    y_at <- pretty(y_all, n = 5)
    axis(1, at = x_at, labels = paste0(round(x_at * 100), "%"))
    axis(2, at = y_at, labels = paste0(round(y_at * 100), "%"))
    }
  
  # block: block rate vs PG
  pg_scatter(
    all_data  = pg_block,
    harv_data = pg_block,
    x_col = "pct_is_block", y_col = "pg",
    xlab  = "Block Rate (min. 30 opportunities)",
    ylab  = "Points Gained per Opportunity",
    main  = "Block",
    min_n_col = "n", min_n = 30
  )
  
  # set: kill rate of sets vs PG
  pg_scatter(
    all_data  = pg_set,
    harv_data = pg_set,
    x_col = "pct_kill_no_block", y_col = "pg",
    xlab  = "Kill Rate of Sets - No Block (min. 100 sets)",
    ylab  = "Points Gained per Opportunity",
    main  = "Set",
    min_n_col = "n", min_n = 100
  )
  
  dev.off()
  cat("Scatter plots saved to ~/Downloads/pg_scatter_harvard.pdf\n")
}
