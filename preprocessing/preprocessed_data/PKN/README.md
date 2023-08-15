## Prior Knowledge Network (PKN)
Files resulting from building a prior knowledge network.

## Files in this folder 
1. add.tsv : Interactions to manually add in the network.
2. edit.tsv : Interactions to manually edit in the network.
3. remove.tsv : Interactions to manually remove from the network.
4. L1000_latest_Add_lvl3.tsv : Some transcription factors (TFs) are directly targeted by drugs and so we will try to forcefully keep them in the network when trimming it (IF possible).
5. L1000_lvl3_DT.tsv : Drug-target interactions in long format. They will be treated also like ligand-receptor interactions.
6. pknFull.tsv : Initial PKN from retrieving interactions from OmniPath and manual curation.
7. pkn.tsv : Initial PKN after filtering interactions where there is disagreement about the sign of the interaction between multiple sources.
8. RLFull.tsv : Extracted receptor-ligand interactions from OmniPath.
9. RL.tsv : Extracted receptor-ligand interactions from OmniPath, after very short filtering and keeping only trusted sources.
19. l1000_lvl3_withsignor-Annotation.tsv : Annotation of the final PKN signaling model, after trimming, which will be used to train machine learning models.
11. l1000_lvl3_withsignor-Model.tsv : Final PKN signaling model, after trimming, which will be used to train machine learning models.
