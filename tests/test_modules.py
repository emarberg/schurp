from modules import QPModule


def test_create():
    m = QPModule.create_gelfand_a(0, 0)
    m = QPModule.create_gelfand_a(1, 0)
    m = QPModule.create_gelfand_a(2, 0)
    m = QPModule.create_gelfand_a(2, 1)
    m = QPModule.create_gelfand_a(3, 1)
