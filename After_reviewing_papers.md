# After-Reading Paper Notes: Points, Thoughts, Questions, etc.

## Adaptive immunity to SARS-CoV-2 and COVID-19
[Adaptive immunity to SARS-CoV-2 and COVID-19](https://pubmed.ncbi.nlm.nih.gov/33497610/)

### Points
- SARS2-specific T-cell responses are significantly associated with the milder disease — T-cell responses are important for control and resolution of infection.
- The inability of the innate immune system to produce effective IFN (interferon) has been closely linked to the inability to control an initial SARS2 infection.
- CD4+ T response to SARS2 > CD8+ T and has been associated with control of primary SARS2 infection.
- CD4+ T cell s had the strongest association with lessened covid19 severity (compared with CD8+ T).
- Different T and B cell reactions might result from a disconnect between B cell and T cell responses, stemming from changes in early innate immune responses. This could potentially lead to delays in kinetics or dysregulation of T cell priming.
- A practical working model suggests that the severity and duration of COVID-19 are primarily determined by the evasion of early innate immune recognition and the subsequent kinetics of the adaptive immune response.
- It’s notable that the an absence of a correlation between neutralizing antibodies and recovery from covid.
- The ability to control COVID without substantial contribution from neutralizing antibodies, as long as a strong T cell response is present.

### Thoughts & Questions
- Why doesn't the innate response decrease in severe SARS2 infections? to compensate for T cells? miss a connection with T-cell response (also mentioned above)
- As mentioned in the paper "Development of neutralizing antibodies against SARS2 is relatively easy, because it can be accomplished by many B cells with little or no affinity maturation required", "SARS2 neutralizing antibody response generally develop from **naive B cells**", "neutralizing epitopes...easy to recognized by antibodies". Then what causes the late antibody generation in severe cases? must have correlations.
- Heterogeneity in adaptive immunity (would be interesting to investigate). From my perspective, I don't think pre-existing immunity would contribute a lot to the heterogeneity, considering myself as an example.

## Dawn Bowdish Grant
[Dawn Bowdish grant](https://mcmasteru365-my.sharepoint.com/:f:/g/personal/xien6_mcmaster_ca/EgAXnNQUnjdHr5uDwCUVkucBKblk85ETDEA_jaEdvbnR2Q)

### Points
- Neutralizing antibodies are not a strong correlate of protection in older adults.
- Non-neutralizing antibodies and cellular immunity are the correlates of protection.
- Antibodies capable of mediating FC-mediated antibody functions such as ADP and ADCC were generated after the first dose of vaccine, possibly contributing to the production.
- Differences in antibody function may explain the protectiveness of different vaccinations (studied in aim 1)
- Hypothesis: non-neutralizing antibody + cellular immune response (caused by different vaccines) + hybrid immune —> contribute to susceptibility and resistance to infection in older adults.
- Increasing the amount of antigen in a vaccine will increase immunogenicity in older adults (??).
- Having a ‘mixed’ vaccine regime is associated with an increased risk of infection (even after 4th and 5th doses) (too general?).
- Hypothesis: Having had a primary series of Modern is associated with lasting differences in cellular and humor immune responses that may explain the prolonged protection (studied in aim3).
- COVID-19 infection did not always lead to protection from subsequent infections
- The infection risk increased in the BA.5 wave in individuals who had had an infection in the first Omicron wave
- Certain SARS2 variants may increase infection risk by altering leukocyte number or function (first found)

### Thoughts & Questions
- How to define **"infection"** in the context of this paper?
- Consider the factor (neutralizing antibodies, GM-CSF expressing Spike-specific CD8+ T cells, etc.) separately, not as an integral or digging their relationships (but we already talked about that as not the main goal here)
- Mentioned the effect of 'mixed vaccine', but different vaccines have different compositions (other than the antigens). Should we consider the interactions between 'other substances' in the vaccine and their impact on the subjects?
- Is there any trade-offs between the protection acquired from multiple vaccinations (or more antigen, as an initial step considering all the vaccine are from the same company) and adverse effects (side effects are not mentioned in the paper)
- Moderna: Adverse effects may potentially contribute to the onset of additional health conditions in the elderly population (not mentioned)  (??)
- The order of experiments 3.1, 3.2


## Immune boosting bridges leaky and polarized vaccination models (our paper)
[Immune boosting bridges leaky and polarized vaccination models](https://www.medrxiv.org/content/10.1101/2023.07.14.23292670v2)
### Points
- **Leaky vaccination model:** assumes all vaccinated individuals experience a reduced force of infection by the same amount
  - Drawback: unrealistic assumption that vaccinated individuals who are exposed to infection can still remain susceptible, independent of previous exposures
- **Polarized vaccination Model:** Some factions of vaccinated individuals are completely protected, while the remaining fraction remains completely susceptible
  - Drawback: Extreme assumption that a fraction of vaccinated individuals do not receive any protection
- Outbreak size (P) <  (L)
- **Immune boosting Model:** Vaccinated, yet susceptible, individuals can gain protection without developing a transmissible infection
- The boasting model predicts identical epidemic dynamics as the polarized vaccination model (build connection)
- **Generalized vaccination model:** Combine the previous three models -- explore how the assumption of immunity affects epidemic dynamics and estimates of vaccine effectiveness
- Different outcomes: for a high force of infection, almost all individuals eventually get infected (L), whereas many individuals are permanently protected in the polarized model.
- Fact that the leaky model doesn’t consider: vaccinated individuals who successfully fight exposures (no infectious or developing symptoms) can experience immune boosting -> become less susceptible to future infections
- Leaky vaccine — history based: susceptibility of an individual depends only on their history of infections (?)
- Polarized vaccine — status based: keep tract of immune statuses of individuals, rather than their infection histories (?)
- Immune boosting model assumes that unsuccessful challenges elicit immune response -> move individuals from S(v) to R(v) at rate VE(L)λ(t)
- Epidemic dynamics are independent of the shape of the *susceptibility distribution* under immune boosting; epidemic dynamics are sensitive to the susceptibility distribution under a leaky vaccination model (?)
- Two ways of estimating vaccine effectiveness: cumulative incidence & instantaneous hazard (?)
- Cumulative incidence based effectiveness estimates will reflect initial efficacy for polarized vaccination and immune boosting models; Hazard based estimates reflect efficacy for the leaky vaccination model

### Thoughts & Questions
- Consider the n^th infection by adding a new parameter
- Why assume a different recovery rate?
- History/status bases?
- Susceptibility distribution and what impact (beneficial?) does this have on the model's predictions?
- Not quite understand why define immune-status as R(u), S(v), R(v) in leaky model
- *logics behind the formula* — Generalized model
- Equations: (26-33) 
