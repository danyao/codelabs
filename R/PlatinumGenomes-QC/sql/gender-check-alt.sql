# Compute the inbreeding coefficient F for each individual within chromosome X to
# help determine whether the gender phenotype value is correct for each individual.
# F >= 0.5 ==> Male
# F <  0.2 ==> Female

SELECT
  individual,O_hom, E_home, N_sites, (O_hom - E_hom) / (N_sites - E_hom) as F
FROM (
  SELECT
    individual,
    COUNT(*) as N_sites,
    SUM(first_allele == second_allele) as O_hom,
    SUM(1.0 - 2.0 * freq * (1.0 - freq) *
        (called_allele_count / (called_allele_count - 1.0))) AS E_hom,
  FROM (
    SELECT
      call.call_set_name as individual,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
      SUM(call.genotype >= 0) WITHIN RECORD AS called_allele_count,
      IF(SUM(call.genotype = 1) > 0,
         SUM(call.genotype = 1) / SUM(call.genotype >= 0),
         -1) WITHIN RECORD AS freq,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts
    FROM [google.com:biggene:platinum_genomes.expanded_variants]
    WHERE
      reference_name = 'chrX'
      AND reference_bases in ('A', 'C', 'G', 'T')
      # Omit pseudoautosomal regions
      AND start NOT BETWEEN 59999 AND 2699519
      AND start NOT BETWEEN 154931042 AND 155260559
    OMIT call IF
      # Skip no calls and haploid sites.
      SOME(call.genotype < 0) OR COUNT(call.genotype) != 2
    HAVING
      # Only consider entries where one alternate base exists.
      num_alts = 1
      # Skip entries where all samples have the same allele.
      AND freq > 0 AND freq < 1
  )
  GROUP BY individual
)
ORDER BY individual
