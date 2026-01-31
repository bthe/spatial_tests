library(tidyverse)

df <-
  readr::read_csv('ldists.csv') |>
  dplyr::filter(sampling_type == 30, year %in% conf$years) |>
  dplyr::select(
    year,
    month,
    station = sample_id,
    begin_lat,
    begin_lon,
    depth,
    tow_length,
    tow_start,
    tow_end,
    length,
    N
  ) |>
  mutate(
    length = pmax(
      pmin(
        conf$maxLength,
        conf$dLength * floor(length / conf$dLength)
      ),
      conf$minLength
    ),
    duration = as.numeric(tow_end - tow_start) / 60,
  ) |>
  group_by(
    year,
    month,
    station,
    begin_lat,
    begin_lon,
    Y = begin_lat,
    X = begin_lon,
    depth,
    tow_length,
    tow_start,
    tow_end,
    length,
    duration
  ) |>
  summarise(N = sum(N, na.rm = TRUE))

attr(df, 'projection') <- "LL"
attr(df, 'zone') <- 27
altMax = -999
altMin = 999

df <-
  df |>
  PBSmapping::convUL() |>
  dplyr::rename(UTMX = X, UTMY = Y) |>
  ungroup() |>
  sun_iceland_mutate(
    tow_start,
    begin_lat,
    begin_lon,
    prefix = ''
  ) |>
  rename(sunAlt = altitude_deg) |>
  mutate(
    sunsetting = tow_start < sunrise_utc_hours | tow_start > sunset_utc_hours,
    # sunsetting = ifelse(sunsetting, 1, 0),
    dist = 4,
    timeInDay = lubridate::hour(tow_start) + lubridate::minute(tow_start) / 60,
    sunAltTrans = ifelse(
      sunsetting,
      (sunAlt - altMin) / (-altMin + altMax) / 2,
      1 - (sunAlt - altMin) / (-altMin + altMax) / 2
    ),
    sunAltTrans = pmin(pmax(sunAltTrans, 0), 1)
  ) |>
  select(
    Station = station,
    LengthGroup = length,
    #LengthInterval = length_interval,
    fishObs = N,
    UTMX,
    UTMY,
    dist,
    timeInDay,
    x = begin_lon,
    y = begin_lat,
    year,
    sunAlt,
    sunAltTrans,
    #sunsetting,
    date = tow_start,
    depth
  )

# df |>
#   split(df$year) -> data
