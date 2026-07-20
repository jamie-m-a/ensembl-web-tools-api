import unittest
from pydantic import ValidationError

from vep.models.vcf_results_model import (
    Annotation,
    PaginationMetadata,
    PredictedIntergenicConsequence,
    PredictedTranscriptConsequence,
    FeatureType,
    Strand,
    ReferenceVariantAllele,
    Location,
    AlternativeVariantAllele,
    Variant,
    Metadata,
    VepResultsResponse,
)


class TestVCFResultModel(unittest.TestCase):

    def test_predicted_intergenic_consequence(self):
        consequence = PredictedIntergenicConsequence()
        self.assertIsNone(consequence.feature_type)
        self.assertEqual(consequence.consequences, ["intergenic_variant"])

    def test_predicted_transcript_consequence(self):
        consequence = PredictedTranscriptConsequence(
            feature_type=FeatureType.transcript,
            stable_id="ENST00000367770.8",
            gene_stable_id="ENSG00000157764.13",
            gene_symbol=None,
            biotype="protein_coding",
            is_canonical=True,
            consequences=["missense_variant"],
            strand=Strand.forward,
        )
        self.assertIsNone(consequence.gene_symbol)

        # Valid instance without explicitly setting gene_symbol (default=None)
        consequence_no_symbol = PredictedTranscriptConsequence(
            feature_type=FeatureType.transcript,
            stable_id="ENST00000367770.8",
            gene_stable_id="ENSG00000157764.13",
            biotype="protein_coding",
            is_canonical=True,
            consequences=["missense_variant"],
            strand=Strand.forward,
        )
        self.assertIsNone(consequence_no_symbol.gene_symbol)

    def test_alternative_variant_allele(self):
        alternative_allele = AlternativeVariantAllele(
            allele_sequence="A",
            allele_type="insertion",
            colocated_variants=["rs123"],
            annotations=[
                Annotation(plugin="spdi", scope="allele", data={"spdi": "1:1:A:T"})
            ],
            predicted_molecular_consequences=[
                PredictedIntergenicConsequence()
            ],
        )
        self.assertEqual(alternative_allele.colocated_variants, ["rs123"])
        self.assertEqual(alternative_allele.annotations[0].plugin, "spdi")

        # Valid instance without the optional annotation fields (both default
        # to empty lists)
        alternative_allele_bare = AlternativeVariantAllele(
            allele_sequence="A",
            allele_type="insertion",
            predicted_molecular_consequences=[
                PredictedIntergenicConsequence()
            ],
        )
        self.assertEqual(alternative_allele_bare.colocated_variants, [])
        self.assertEqual(alternative_allele_bare.annotations, [])

    def test_variant(self):
        variant = Variant(
            name=None,
            allele_type="SNP",
            location=Location(region_name="1", start=10000, end=10001),
            reference_allele=ReferenceVariantAllele(allele_sequence="C"),
            alternative_alleles=[
                AlternativeVariantAllele(
                    allele_sequence="T",
                    allele_type="SNP",
                    predicted_molecular_consequences=[
                        PredictedIntergenicConsequence()
                    ]
                )
            ]
        )
        self.assertIsNone(variant.name)

        # Valid instance without explicitly setting name (default=None)
        variant_no_name = Variant(
            allele_type="SNP",
            location=Location(region_name="1", start=10000, end=10001),
            reference_allele=ReferenceVariantAllele(allele_sequence="C"),
            alternative_alleles=[
                AlternativeVariantAllele(
                    allele_sequence="T",
                    allele_type="SNP",
                    predicted_molecular_consequences=[
                        PredictedIntergenicConsequence()
                    ]
                )
            ]
        )
        self.assertIsNone(variant_no_name.name)

    def test_metadata_required(self):
        metadata = Metadata(
            pagination=PaginationMetadata(page=1, per_page=10, total=100)
        )
        self.assertEqual(metadata.pagination.page, 1)

        # Invalid instance without pagination - should raise ValidationError
        with self.assertRaises(ValidationError):
            Metadata()

    def test_vep_results_response(self):
        metadata = Metadata(
            pagination=PaginationMetadata(page=1, per_page=10, total=100)
        )
        variant = Variant(
            name=None,
            allele_type="SNV",
            location=Location(region_name="1", start=10000, end=10001),
            reference_allele=ReferenceVariantAllele(allele_sequence="A"),
            alternative_alleles=[
                AlternativeVariantAllele(
                    allele_sequence="T",
                    allele_type="SNP",
                    predicted_molecular_consequences=[
                        PredictedIntergenicConsequence()
                    ]
                )
            ],
        )
        response = VepResultsResponse(metadata=metadata, variants=[variant])
        self.assertEqual(response.metadata.pagination.page, 1)
        self.assertEqual(len(response.variants), 1)

        # Invalid instance without metadata - should raise ValidationError
        with self.assertRaises(ValidationError):
            VepResultsResponse(variants=[variant])
