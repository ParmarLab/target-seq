the_plan <-
  drake_plan(
    f1 = figure1(outdir = "plots/", hpsc.file = "processed_data/hpsc.human_harmony.rds",
                 fetal.file = "processed_data/fetal.human_harmony.rds" ),
    f2 = figure2(outdir = "plots/", hpsc.neuron.file = "processed_data/hpsc.human_harmony.neurons.rds",
                 fetal.file = "processed_data/fetal.human_harmony.rds",
                 hpsc.thneuron.file = "processed_data/hpsc.human_harmony.THneurons.rds"),
    f3 = figure3(outdir = "plots/", hpsc.thneuron.file="processed_data/hpsc.human_harmony.THneurons.rds"),
    f4 = figure4(outdir = "plots/", hpsc.thneuron.file = "processed_data/hpsc.human_harmony.THneurons.rds"),
    f5 = figure5(outdir = "plots/", hpsc.thneuron.file="processed_data/hpsc.human_harmony.THneurons.rds")
  )
  