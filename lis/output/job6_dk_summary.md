# job 1506296 gb_job6_dk_floor_decomp (Denmark, one-person households), received 2026-09-27 18:04Z, 146 s

- pubonly (hipubsoc>=0.9 hitotal) among 65+: 0.1-0.4% -> public pension is NOT in hipubsoc; decomposition must use pension components (hi3x).
- old_all w5000 / w10000: 
  04 p=1 ; 07 p=.87 (valley 92.5-102.5, m1 95k m2 130k); 10 p=.09/1 (112.5-132.5, m1 110k, m2 150k); 13 p=.105/1 (127.5-137.5, m1 125k m2 155k);
  16 p=.01/.02 (132.5-147.5 | 125-145, m1 135k/130k, m2 170k); 19 p=.01/.005 (137.5-152.5 | 135-155, m1 140k, m2 180k)
  m1/FP gross: .83 .88 .90 .88 (2010-2019). pslo .10-.15, pshi .23-.27
- u30 trim .01: 10,13,16,19 reject via near-zero spike (m1=0). trim .05/.10: 04 p=.135/.955; 07 .94/.85; 10 .94/.87; 13 .235/.225; 16 .005/.005 (72.5-97.5, m1 65k, m2 100k); 19 .025/.015 (82.5-102.5, m1 65k, m2 105k)
- mid 30-64: trim .01 rejects 16,19 via zero spike; trim .05: 04 1, 07 .935, 10 .975, 13 .72, 16 .06, 19 .145; trim .10: 1,.93,.95,.68,.32,.14
- SU placebo (gross, u30) excess at 0.8/0.9/1/1.1/1.2 x SU:
  04 -.217 -.389 .666 -.106 .079 ; 07 -.096 -.264 .399 .24 -.093 ; 10 -.137 -.436 .558 .339 .022 ; 13 -.281 -.521 .533 .16 -.086 ; 16 -.246 -.353 .747 .083 -.065 ; 19 -.183 -.351 .712 .224 -.141
- FP gross bunching (old) at .8/1/1.2: none systematic (04 -1.49 -.57 -.24; 16 .07 -.09 .01; 19 -.04 .07 .02)
- pubonly dhi quantiles show repeated tiny values (1180, 1230) -> bottom-coded/near-zero records; justifies 5% trim.
