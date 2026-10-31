library(tidyverse)
library(lme4)
library(emmeans)

# =================================================================
# 0. Quali componenti combinare
# =================================================================
# Basta cambiare questo vettore per usare più (o meno) PC, es. 1:3, 1, ma anche c(1, 3)
pcs        <- 1:2
score_vars <- paste0("s", pcs)   # nomi delle colonne degli score in PCscores
                                  # -> assumo "s1", "s2"; correggi se diverso

crit_val <- qnorm(0.975)         # ~1.96, usato per ricostruire la banda di
                                  # confidenza combinata (approssimazione
                                  # normale — vedi nota sotto)

# =================================================================
# 1. Un modello misto per ciascuna PC, poi emmeans in formato "long"
# =================================================================
# Stessa formula del tuo modello originale, solo la variabile risposta cambia
fit_one_pc <- function(score_var) {
  form <- as.formula(
    paste(score_var, "~ dialect * syllable + (1 | word) + (syllable | speaker)")
  )
  lmer(form, data = PCscores)
}

pc_models <- set_names(score_vars, pcs) %>%
  map(fit_one_pc)

# per ogni PC estraggo emmean, SE, lower/upper CL; PC = indice della componente
pc_emmeans <- pc_models %>%
  imap_dfr(~ emmeans(.x, ~ syllable | dialect) %>%
              as_tibble() %>%
              mutate(PC = as.integer(.y)))

# Se ti servono anche i contrasti a coppie per ciascuna PC (come nel tuo
# `pairwise ~ syllable|dialect`), puoi estrarli allo stesso modo per singolo
# modello, es.: emmeans(pc_models$`1`, pairwise ~ syllable|dialect)$contrasts

# =================================================================
# 2. Ricostruzione delle curve: mean + sum_k( emmean_k * eigenfunction_k )
# =================================================================
predCurves <- pc_emmeans %>%
  expand_grid(Dim = DimCurves) %>%
  group_by(syllable, dialect, Dim) %>%
  reframe({

    # IMPORTANTE: con più PC il gruppo ha più di una riga, quindi "Dim"
    # (colonna di raggruppamento) qui dentro sarebbe un vettore tipo c(1,1)
    # invece che uno scalare. mfpca$meanFunction[[c(1,1)]] farebbe
    # un'indicizzazione RICORSIVA (equivalente a mfpca$meanFunction[[1]][[1]])
    # anziché selezionare l'elemento 1 -> lo forzo a scalare qui:
    Dim <- Dim[1]

    # griglia temporale + funzione media, comune a tutte le PC di questo Dim
    base <- funData2long1(mfpca$meanFunction[[Dim]], value = "y0")

    # valori dell'autofunzione della PC "pc"
    eig_y <- function(pc) funData2long1(mfpca$functions[[Dim]][pc], value = "y")$y

    # contributo alla curva media: somma pesata delle autofunzioni
    contrib <- map2(PC, emmean, function(pc, m) {
      m * eig_y(pc)
    }) %>% reduce(`+`)

    # contributo alla varianza punto-per-punto (assumendo PC indipendenti):
    # Var(t) = sum_k( SE_k^2 * eigenfunction_k(t)^2 )
    var_contrib <- map2(PC, SE, function(pc, se) {
      (se * eig_y(pc))^2
    }) %>% reduce(`+`)

    sd_curve <- sqrt(var_contrib)

    base %>%
      mutate(
        y  = y0 + contrib,
        yl = y - crit_val * sd_curve,
        yu = y + crit_val * sd_curve
      ) %>%
      select(-y0)
  })

# =================================================================
# 3. Ricostruzione delle durate (lograte combinato dalle stesse PC)
# =================================================================
predDur <- pc_emmeans %>%
  expand_grid(landmarks2long(reg$landmarks)) %>%
  group_by(syllable, dialect, leftBoundary) %>%
  reframe(
    duration = lograte2duration(
      lograte = mfpca$meanFunction[[DimLograte]] +
        reduce(map2(PC, emmean, ~ .y * mfpca$functions[[DimLograte]][.x]), `+`),
      from = from[1], to = to[1]
    )
  )

# aggiungo la posizione cumulativa dei landmark (come nello script originale)
predDur <- predDur %>%
  group_by(syllable, dialect) %>%
  mutate(land = cumsum(duration))

# =================================================================
# 4. Da qui in poi: IDENTICO al tuo script originale
#    (allineamento temporale, landmark, plot — non dipende dal numero di PC)
# =================================================================
predCurves <- predCurves %>%
  group_by(syllable, dialect) %>%
  mutate(origTime = {
    landmarks <- c(0,
                   predDur %>%
                     inner_join(cur_group(), by = c("dialect", "syllable")) %>%
                     pull(land))
    landmarkreg_timeSamples(timeSamples = time,
                             inputMarks = reg$landmarks,
                             targetMarks = landmarks)
  })

predLand <- predCurves %>%
  ungroup() %>%
  select(dialect, syllable, Dim, origTime, y) %>%
  nest(.by = c(dialect, syllable, Dim), .key = "curves") %>%
  inner_join(
    predDur %>%
      group_by(dialect, syllable) %>%
      summarise(land = list(c(0, land)), .groups = "drop"),
    by = c("dialect", "syllable")
  ) %>%
  mutate(y = map2(curves, land,
                  ~ approx(.x$origTime, .x$y, xout = .y, rule = 2)$y)) %>%
  select(-curves) %>%
  unnest(cols = c(land, y))

# check: should give 0 rows, otherwise something went wrong

predCurves %>% ungroup() %>% count(dialect, syllable, Dim, origTime) %>% filter(n > 1)

# change DIM names
predCurves = predCurves %>% rename(Channel = Dim)
# Rename levels
predCurves$Channel <- plyr::revalue(as.character(predCurves$Channel), c("1" = "oral", "2" = "nasal"))

#plot everything onto one panel

#syll.colors <- c("purple", "slategrey")

DimLabels <- c(`1` = "oral", `2` = "nasal")

ggplot(predCurves) +
  aes(origTime, y, linetype = syllable, color = dialect) +
  geom_line(linewidth = 1.2) +
  # geom_point(data = predLand,
  #            mapping = aes(land, y, color = dialect),
  #            inherit.aes = FALSE,
  #            size = 3
  # ) +
  facet_grid(~Channel) +
  scale_color_manual(values = c("gold3", "forestgreen"),labels = c("NMF", "QF")) +
  #scale_color_manual(values = syll.colors) +
  xlab("Time (seconds)") +
  ylab("Pa") + 
  mytheme +
  theme(legend.position = "top",
        strip.text.x = element_text(color = "black", size = 12),
        strip.text.y = element_text(color = "black", size = 12))

# =================================================================
# 5. New part compared to ICPhS analysis: calculate nasalance
# from *reconstructed* oral and nasal curves
# =================================================================
#
# save df before transforming it

write.table(predCurves, "predCurves_oral-nasal.txt")
# predCurves <- read.table("predCurves_oral-nasal.txt", header = TRUE)

dplyr::count(predCurves, syllable, dialect, Channel)   # 100 time points for each combination

# 2. Pivot long to wide (oral e nasal side by side) --------
# one row = one syllable x dialect x time point - combination
nas <- predCurves %>%
  select(syllable, dialect, Channel, time, origTime, y) %>%
  pivot_wider(names_from  = Channel,
              values_from = y) %>%
  mutate( nasalance = nasal / (nasal + oral) # calculates nasalance by each timepoint
  )

#save for the future
write.csv(nas, "nasalance_curves.csv", row.names = FALSE)

# Controllo: nessun NA (cioe' oral e nasal allineati su ogni time point)
stopifnot(!anyNA(nas$nasalance))
head(nas)

# Plot: all 4 curves in a unique panel
p_all <- ggplot(nas,
  aes(origTime, nasalance, linetype = syllable, color = dialect)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("gold3", "forestgreen"),labels = c("NMF", "QF")) +
   # scale_y_continuous(labels = scales::percent)  #if you want nas in percentages
  xlab("Time (seconds)") +
  ylab("Nasalance") + 
  mytheme +
  theme(legend.position = "top",
  strip.text.x = element_text(color = "black", size = 12),
  strip.text.y = element_text(color = "black", size = 12))+
  ggtitle ("Nasalance curves predicted by MFPCA on oral and nasal channels")
                
p_all

# 6. Salvataggio ---------------------------------------------

ggsave("nasalance.png", p_all, width = 9, height = 4.5, dpi = 300) # change sets as you wish


