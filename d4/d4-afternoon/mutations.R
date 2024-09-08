library(tidyverse)
library(broom)

dnm <- read_csv(file = "~/qbb2024-answers/d4/d4-afternoon/aau1043_dnm.csv")

ages <- read_csv(file = "qbb2024-answers/d4/d4-afternoon/aau1043_parental_age.csv")

dnm_summary <- dnm %>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm = sum(Phase_combined == "father", na.rm = TRUE),
            n_maternal_dnm = sum(Phase_combined == "mother", na.rm = TRUE))

joined_list <- left_join(dnm_summary, ages, by = "Proband_id")

ggplot(data = joined_list,
        mapping = aes(x = Mother_age , y = n_maternal_dnm)) +
  geom_point()

ggplot(data = joined_list,
       mapping = aes(x = Father_age , y = n_paternal_dnm)) +
  geom_point()

mother_model <- lm(data = joined_list,
   formula = n_maternal_dnm ~ 1 + Mother_age)
  summary()

father_model <- lm(data = joined_list,
   formula = n_paternal_dnm ~ 1 + Father_age)
  summary()
#There is a significant positive relationship between the two variables in BOTH the maternal and paternal cases which matches what I would have predicted based on the plot.
#The p-value for both is well below 0.05, indicating that both relationships are statistically very significant.
#The maternal slope and R-squared indicate that, while there is a relationship between the variables, it is neither very large (0.38) nor strong (0.23).
#The paternal slope and R-squared, on the other hand, show a very large (1.35) and strong (0.62) relationship.
 
new_dude <- tibble(Father_age = 50.5)
predict(father_model, newdata = new_dude)
#By generating a new father based on the model above we can predict that the number of paternally inherited dnms will be around 78.

#Question 2.5
ggplot() +
  geom_histogram(data = joined_list, mapping = aes(x = n_maternal_dnm, fill = "Maternal"), alpha = 0.5) +
  geom_histogram(data = joined_list, mapping = aes(x = n_paternal_dnm, fill = "Paternal"), alpha = 0.5) +
  labs(x = "Maternal and Paternal DNMs", fill = NULL)

#Question 2.6
n_maternal_dnm <- joined_list$n_paternal_dnm
n_paternal_dnm <- joined_list$n_maternal_dnm
t.test(n_maternal_dnm, n_paternal_dnm, var.equal = FALSE)
#A Welch's t-test is appropriate here because we are comparing one continuous variable across two independent groups.
#The test is also optimal because of the difference in variance between the groups, their normal distribution, and our large sample size.
#The p-value indicates that the difference between the groups is very significant.
#The t-value is also high, indicating a large difference between the means of the two groups.
#Altogether the results show that human offspring are significantly more likely to inherit DNMs paternally. 






