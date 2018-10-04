# Configuration de test

Processeur de portable i7-4720HQ (4 coeurs + HT)
Rust 1.29.1
Commit cc1f0481d6660e6aa5581f68f2f0416e88d09bec
Running on 10^7 events.


# Tests de features isolées

Mode par défaut :
- Exactement reproductible
- 5s -> 3.4s (32% plus rapide qu'original)

Mode faster-evgen :
- Différences de ~0.1% vs défaut (sauf résultat très instable I_MX)
- 2.7s (21% vs défaut, 46% vs original)

Mode f32 :
- Différences de ~0.1% vs défaut
- 2.5s (26% vs défaut, 50% vs original)

Mode multi-threading :
- Exactement reproductible
- 0.81s (speedup 4.2x vs défaut)

Modes multi-threading + faster-threading :
- Stochastique + différences de ~0.1%
- 0.79s (différence négligeable avec mt seul)

Mode no-photon-sorting :
- Exactement reproductible
- 3.4s (= défaut)

Mode standard-random :
- Différences de 0.1% vs défaut
- 3.4s (= défaut)


# Tests de features combinées

Modes faster-evgen + f32 :
- Différences de 0.1% vs défaut
- 2.11s (40% vs défaut, 67% vs original)

Modes faster-evgen + multi-threading :
- Exactement reproductible vs faster
- 0.87s (speedup 3.1x => Bottleneck !)

Modes faster-evgen + multi-threading + faster-threading :
- Différences de 0.1% vs faster
- 0.57s (speedup 4.7x => OK)

Modes faster-evgen + multi-threading + faster-threading + standard-random :
- Différences de 0.1% vs faster
- 0.59s (speedup 4.6x => OK)

Modes faster-evgen + f32 + multi-threading + faster-threading :
- Différences de 0.1% vs défaut
- 0.48s (speedup 4.4x => OK)