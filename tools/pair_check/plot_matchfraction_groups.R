library("readxl")
library("tidyverse")


match <- as.data.frame(read_excel("ws_parcheck_match.xlsx"))
no_match <- as.data.frame(read_excel("ws_parcheck_nonmatch.xlsx"))
family <- as.data.frame(read_excel("ws_parcheck_family.xlsx"))


match$`Match Fraction` <- as.numeric(match$`Match Fraction`)
no_match$`Match Fraction` <- as.numeric(no_match$`Match Fraction`)
family$`Match Fraction` <- as.numeric(family$`Match Fraction`)

match <- match %>% add_column(Group = 'Match')
no_match <- no_match %>% add_column(Group = 'Mismatch')
family <- family %>% add_column(Group = 'Family')

match <- match %>% add_column(dnorm = dnorm(match$`Match Fraction`, mean = mean(match$`Match Fraction`), sd = sd(match$`Match Fraction`)))
no_match <- no_match %>% add_column(dnorm = dnorm(no_match$`Match Fraction`, mean = mean(no_match$`Match Fraction`), sd = sd(no_match$`Match Fraction`)))
family <- family %>% add_column(dnorm = dnorm(family$`Match Fraction`, mean = mean(family$`Match Fraction`), sd = sd(family$`Match Fraction`)))





match_fractions <-
  full_join(match %>% dplyr::select(`Match Fraction`, Group, dnorm),
            no_match %>% dplyr::select(`Match Fraction`, Group, dnorm))
match_fractions <- full_join(match_fractions, family %>% dplyr::select(`Match Fraction`, Group, dnorm))



pl_all <- ggplot(match_fractions, aes(`Match Fraction`, dnorm)) +
  geom_point(aes(color = Group)) +
  scale_x_continuous(breaks = seq(0,1,by=0.05)) + 
  ylab("Probability density") +
  geom_vline(xintercept = mean(match$`Match Fraction`)+sd(match$`Match Fraction`), linetype="dotted", color = "green") + # one sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)-sd(match$`Match Fraction`), linetype="dotted", color = "green") + # one sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)+2*sd(match$`Match Fraction`), linetype="dotted", color = "green") + # two sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)-2*sd(match$`Match Fraction`), linetype="dotted", color = "green") + # two sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)+3*sd(match$`Match Fraction`), linetype="dotted", color = "green") + # three sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)-3*sd(match$`Match Fraction`), linetype="dotted", color = "green") + # three sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)+sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # one sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)-sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # one sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)+2*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # two sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)-2*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # two sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)+3*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # three sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)-3*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # three sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)+sd(family$`Match Fraction`), linetype="dotted", color = "red") + # one sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)-sd(family$`Match Fraction`), linetype="dotted", color = "red") + # one sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)+2*sd(family$`Match Fraction`), linetype="dotted", color = "red") + # two sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)-2*sd(family$`Match Fraction`), linetype="dotted", color = "red") + # two sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)+3*sd(family$`Match Fraction`), linetype="dotted", color = "red") + # three sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)-3*sd(family$`Match Fraction`), linetype="dotted", color = "red") # three sd from mean


pdf("groups_match_fraction.pdf", width = 15, height = 7.5)
plot(pl_all)
dev.off()

pl_match <- ggplot(match, aes(`Match Fraction`, dnorm)) +
  geom_point(color = "green4") +
  ylab("Probability density") +
  geom_vline(xintercept = mean(match$`Match Fraction`)+sd(match$`Match Fraction`), linetype="dotted", color = "green4") + # one sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)-sd(match$`Match Fraction`), linetype="dotted", color = "green4") + # one sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)+2*sd(match$`Match Fraction`), linetype="dotted", color = "green4") + # two sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)-2*sd(match$`Match Fraction`), linetype="dotted", color = "green4") + # two sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)+3*sd(match$`Match Fraction`), linetype="dotted", color = "green4") + # three sd from mean
  geom_vline(xintercept = mean(match$`Match Fraction`)-3*sd(match$`Match Fraction`), linetype="dotted", color = "green4")  # three sd from mean

pdf("match_plot.pdf", width = 15, height = 7.5)
plot(pl_match)
dev.off()


pl_mismatch <- ggplot(no_match, aes(`Match Fraction`, dnorm)) +
  geom_point(color = "blue") +
  ylab("Probability density") +
  geom_vline(xintercept = mean(no_match$`Match Fraction`)+sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # one sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)-sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # one sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)+2*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # two sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)-2*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # two sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)+3*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue") + # three sd from mean
  geom_vline(xintercept = mean(no_match$`Match Fraction`)-3*sd(no_match$`Match Fraction`), linetype="dotted", color = "blue")  # three sd from mean

pdf("mismatch_plot.pdf", width = 15, height = 7.5)
plot(pl_mismatch)
dev.off()


pl_family <- ggplot(family, aes(`Match Fraction`, dnorm)) +
  geom_point(color = "red") +
  ylab("Probability density") +
  geom_vline(xintercept = mean(family$`Match Fraction`)+sd(family$`Match Fraction`), linetype="dotted", color = "red") + # one sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)-sd(family$`Match Fraction`), linetype="dotted", color = "red") + # one sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)+2*sd(family$`Match Fraction`), linetype="dotted", color = "red") + # two sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)-2*sd(family$`Match Fraction`), linetype="dotted", color = "red") + # two sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)+3*sd(family$`Match Fraction`), linetype="dotted", color = "red") + # three sd from mean
  geom_vline(xintercept = mean(family$`Match Fraction`)-3*sd(family$`Match Fraction`), linetype="dotted", color = "red")  # three sd from mean



pdf("family_plot.pdf", width = 15, height = 7.5)
plot(pl_family)
dev.off()

