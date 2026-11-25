import pytest
from benchmarks.audit_ahif_foster_original import resolve,parse_eml,schema_csv

def test_resolve_exact_version():
    b="https://pasta.edirepository.org/package/"
    a=b+"metadata/eml/knb-lter-kbs/1/2\n"+b+"data/eml/knb-lter-kbs/1/2/"+"a"*32
    assert resolve(a)==("knb-lter-kbs/1/2",["a"*32])
    bad=b+"metadata/eml/knb-lter-kbs/1/2\n"+b+"data/eml/knb-lter-kbs/1/3/"+"a"*32
    assert bad != a
    with pytest.raises(ValueError,match="version"):
        resolve(bad)
def test_eml_and_schema():
    e=b'<eml packageId="edi.1.2"><dataset><title>Hay meadow restoration</title><dataTable><entityName>Biomass</entityName><attributeList><attribute><attributeName>Year</attributeName><attributeDefinition>Harvest year</attributeDefinition></attribute></attributeList></dataTable></dataset></eml>'
    assert parse_eml(e,"edi/1/2")["entities"][0]["attributes"][0]["name"]=="Year"
    z=schema_csv(b"Year,Plot,Seed,Mass\n2000,1,yes,10\n2001,1,yes,12\n")
    assert z["n_rows"]==2 and z["nonrectangular_rows"]==0
