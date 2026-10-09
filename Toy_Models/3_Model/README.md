# 3_Model: Unified 2D Beam/Truss Railway Bridge FEM Model for SHM

Questo modulo unifica e supera le limitazioni dei modelli precedenti, combinando:
1. **La cinematica 2D generale del Modello 1**: travi nel piano con deformazione assiale ($EA$), flessione di Euler-Bernoulli ($EJ$) e rotazione trigonometrica $\mathbf{\Lambda}(\gamma)$.
2. **Il motore dinamico e la modellazione SHM del Modello 2**:
   - Elementi concentrati (molle a terra, masse concentrate, smorzatori).
   - Generazione delle forze da convoglio ferroviario con carichi viaggianti multipli ($n_P$ assi).
   - Analisi modale avanzata con smorzamento e parametri modali.
   - Integrazione temporale diretta ad alte prestazioni con l'algoritmo di **Newmark accelerazione media** ($\beta = 1/4, \gamma = 1/2$) con fattorizzazione di Cholesky precalcolata.
3. **Flessibilità e scenari di danno**:
   - Lettura e parametrizzazione trasparente tramite file `.inp` (`loadstructure.m`).
   - Possibilità di modellare indifferentemente sia **travate da ponte orizzontali** (`bridge_beam.inp`), sia **travature reticolari** (`bridge_truss.inp`), sia archi o pile.
   - Introduzione immediata di scenari di danno (flessionale o assiale) semplicemente modificando l'indice di proprietà dell'asta desiderata nel file `.inp`.

---

## Struttura dei File

| File | Descrizione |
| :--- | :--- |
| `Main.m` | Orchestratore principale della simulazione (caricamento, plot struttura, analisi modale, transito treno, plot accelerazioni). |
| `loadstructure.m` | Parser del file `.inp` (nodi, vincoli, proprietà, travi, molle/masse concentrate e definizione della via di corsa `*DECK`). |
| `assem.m` | Assemblaggio delle matrici globali $M, K, C$ integrando gli elementi continui e i contributi concentrati. |
| `el_tra.m` | Matrici locali a 6 GDL per elemento trave 2D (assiale + flessione) e rotazione al sistema globale. |
| `modal_analysis.m` | Solutore agli autovalori, normalizzazione forme modali e calcolo parametri modali. |
| `Newmark_accelerazioni.m` | Risolutore dinamico passo-passo per il calcolo delle serie temporali di accelerazione e spostamento. |
| `Forze_carichi_viaggianti.m` | Calcolo delle forze nodali equivalenti da convoglio ferroviario su elementi 2D del deck. |
| `dis_stru.m` | Visualizzazione della geometria indeformata, nodi, vincoli e proprietà delle travi. |
| `modes_plot_conc_element.m` | Plot dei profili verticali delle forme modali. |
| `diseg2.m` | Plot della deformata 2D interpolata per le forme modali. |
| `bridge_beam.inp` | File di esempio: travata da ponte a 20 elementi. |
| `bridge_truss.inp` | File di esempio: ponte reticolare Warren con transito sul corrente inferiore. |
