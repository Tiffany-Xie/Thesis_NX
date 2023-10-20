# Proposal

Due 27 October

# Proposed Title
**Within-host dynamical models of acute respiratory viruses: A review of approaches and methods**

# Readings
[Adaptive immunity to SARS-CoV-2 and COVID-19](https://pubmed.ncbi.nlm.nih.gov/33497610/)

[Dawn Bowdish grant](https://mcmasteru365-my.sharepoint.com/:b:/g/personal/xien6_mcmaster_ca/EWOYlawx391Gowrq2UV190UB3B_kw_RTMskuOeI9vxClmg?e=3hXkjt)

[Immune boosting bridges leaky and polarized vaccination models (our paper)](https://www.medrxiv.org/content/10.1101/2023.07.14.23292670v2)

**Linear Chain Trick**

[Olga's thesis](https://macsphere.mcmaster.ca/bitstream/11375/11231/1/fulltext.pdf)

[Olga's thesis (view without download)](https://mcmasteru365-my.sharepoint.com/:b:/g/personal/xien6_mcmaster_ca/Efmc21R_qWROkq4QHVNRNZQBUZsCMQQqdSqahHd9SJr3NQ?e=e9bQ1r)

[Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models](https://link.springer.com/article/10.1007/s00285-019-01412-w)

[Building mean field ODE models using the generalized linear chain trick & Markov chain theory](https://www.tandfonline.com/doi/full/10.1080/17513758.2021.1912418)

# Ideas
- Tropism: What do we know about how or whether tetrapods (or mammals?) have evolved various mechanisms to induce a tradeoff between viral success in lower or upper airways
- Immunity: How do we measure or define immunity? How careful are people in talking about it? We can have protection against infection, the ability to transmit, clinical disease, hospitalization, and death. How do these things wane?

**Mathy stuff:** Simulate, analyze, and fit very simple boxcar models.

# To do (Ningrui)

Let me know what you think of the SARS paper (maybe write a notes document about it on the repo and send me a link or just let me know the name)
* For now, just spend one hour consolidating your _current_ thoughts and questions ✓

Finally, let JD know thoughts about the Dawn Bowdish grant
* ~/Downloads/vaccResponseGrant.pdf ✓
* For now, just spend one hour consolidating your _current_ thoughts and questions ✓

And similarly for [our paper](https://www.medrxiv.org/content/10.1101/2023.07.14.23292670v2)
* For now, just spend one hour consolidating your _current_ thoughts and questions ✓
* Also, work on the math and/or simulation side to see what you can learn ✓

* Set up a oneDrive folder, share it with JD, and add the grant pdf ✓
* This one is very hard, but try to come up with some questions or thoughts ✓

Oct 3, 2023
* Test Erlang using Gamma ✓
* Match the Estimated P to the exact P (approximate using small interval OR *other method*) ✓
* Find better ways to dynamically change the *dI* equations in the model ✓

## 2023 Oct 17 (Tue)

We want to be able to fix the number of boxes and change r and see if that's a better overall approach than changing the number of boxes.

Our starting point is going to be Erlang distributions with n=2, 4 and 8. What values of r match these distributions if we use the pseudo-Erlang with 12 boxes?

We can look both at visual matches and also at shape parameters (for now, we will always keep the mean duration constant).

Our main shape parameter is κ = σ²/μ. We can calculate that analytically from the formulas on the paper. It might be fun to also calculate it numerically from the model output, but that is probably not necessary unless we run into trouble.

After we do visual comparisons, we can also try to go backward. Given a number of boxes and a desired value of κ, it should be possible to numerically find the value of r that matches κ. The best way to do this is probably with the r function `uniroot`. If we do this for κ=1/n_erlang, we can compare pseudo-Erlang and Erlang distributions with the same D and the same κ (we will typically adjust a as the last step to make the values of D match).

2023 Oct 20 (Fri)
=================

For proposal: review a little bit about SInR modeling (linear chains). You can rely mostly on the friendly parts of Olga's thesis. Briefly explain your own understanding of what's going on.

Look to see if Gholami and/or Heffernan has done anything with fitting multiple models with different numbers of boxes. See for example https://pubmed.ncbi.nlm.nih.gov/36773843/, there may also be others on archives?

Find out about our lab's SIR fitting stuff? https://github.com/bbolker/fitsir/

Here is some background about how different kinds of models come together. Do NOT worry about the math details; if the overview helps you, or you find something to cite, that's great. Otherwise you can ignore the whole thing.

Biological parameters are D and κ

Erlang mechanistic parameters are γ and n [clumsy and discrete]

PE mechanistic parameters [fix n]: a and r
