get_fixes <- function(version) {
  if(version == "latest") {
    fixes <- readr::read_csv("hu12_fixes/hu_fixes.csv") %>%
      bind_rows(list(HUC12 = "180102040904", TOHUC = "180102041003", comment = "misdirected"))
  } else if(version == "nhdplusv2") {
  ##### Fixes #####
  fixes <- tibble(HUC12 = character(0),
                  TOHUC = character(0),
                  comment = character(0)) %>%
    bind_rows(list(HUC12 = "101800100703", TOHUC = "101800100702",
                   comment = "self directed")) %>%
    bind_rows(list(HUC12 = "160300020305", TOHUC = "160300020306",
                   comment = "self directed")) %>%
    bind_rows(list(HUC12 = "170401050305", TOHUC = "170401050306",
                   comment = "circular pair")) %>%
    bind_rows(list(HUC12 = "101800110901", TOHUC = "101800110902",
                   comment = "self directed")) %>%
    bind_rows(list(HUC12 = "101800120602", TOHUC = "101800120603",
                   comment = "self directed")) %>%
    bind_rows(list(HUC12 = "100200070304", TOHUC = "100200070307",
                   comment = "self directed")) %>%
    bind_rows(list(HUC12 = "060101051302", TOHUC = "060101051403",
                   comment = "self directed")) %>%
    bind_rows(list(HUC12 = "170402030303", TOHUC = "170402030305",
                   comment = "circular pair")) %>%
    bind_rows(list(HUC12 = "100302030902", TOHUC = "100302030901",
                   comment = "circular pair"))  %>%
    bind_rows(list(HUC12 = "180102040902", TOHUC = "180102041201",
                   comment = "wrong major path")) %>%
    bind_rows(list(HUC12 = "101900180706", TOHUC = "102001010305",
                   comment = "Major route error")) %>%
    bind_rows(list(HUC12 = "101101010705", TOHUC = "101101010707",
                   comment = "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101010707", TOHUC = "101101011003",
                   comment = "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101011003", TOHUC = "101101011004",
                   comment = "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101010906", TOHUC = "101101011004",
                   comment = "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101010809", TOHUC = "101101010707",
                   comment = "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101011002", TOHUC = "101101011003",
                   comment = "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101102050609", TOHUC = "101101013001",
                   comment =  "fix major path")) %>%
    bind_rows(list(HUC12 = "101301010702", TOHUC = "101301010705",
                   comment =  "fix major path")) %>%
    bind_rows(list(HUC12 = "102400110108", TOHUC = "102400110109",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "102400110109", TOHUC = "102400110301",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160302", TOHUC = "170102160308",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160305", TOHUC = "170102160308",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160306", TOHUC = "170102160308",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160308", TOHUC = "170102160000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160000", TOHUC = "170200010203",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170200010203", TOHUC = "170200010205",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900030604", TOHUC = "101900030605",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900030605", TOHUC = "101900030606",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900030505", TOHUC = "101900030606",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900040403", TOHUC = "101900040404",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900040404", TOHUC = "101900030402",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900030406", TOHUC = "101900030407",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900050602", TOHUC = "101900050603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051100010109", TOHUC = "051100010112",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051100010112", TOHUC = "051100010115",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051100010115", TOHUC = "051100010116",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051100010110", TOHUC = "051100010112",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051100010111", TOHUC = "051100010112",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051100010114", TOHUC = "051100010115",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101011003", TOHUC = "101101011004",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101010906", TOHUC = "101101011004",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101011002", TOHUC = "101101011004",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101010707", TOHUC = "101101011004",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101010705", TOHUC = "101101010707",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101101010809", TOHUC = "101101010707",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111401010405", TOHUC = "111401010407",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111401010407", TOHUC = "111401010801",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111401020209", TOHUC = "111401010801",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111401020207", TOHUC = "111401020209",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111401020208", TOHUC = "111401020209",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111401010406", TOHUC = "111401010407",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051401010605", TOHUC = "051401010904",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051401010904", TOHUC = "051401040102",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051401010905", TOHUC = "051401010904",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051401010903", TOHUC = "051401010904",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051401010902", TOHUC = "051401010904",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051401010806", TOHUC = "051401010904",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "140401061008", TOHUC = "140600010303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "140600010303", TOHUC = "140600010306",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "140600010104", TOHUC = "140600010303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "140500020616", TOHUC = "140600010303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "140600010301", TOHUC = "140600010306",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "140600010302", TOHUC = "140600010306",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150502020304", TOHUC = "150502020401",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150502020102", TOHUC = "150502020304",
                   comment =  "UNKNOWN fix (boader hack)")) %>%
    bind_rows(list(HUC12 = "150503010104", TOHUC = "150503010302",
                   comment =  "UNKNOWN fix (boader hack)")) %>%
    bind_rows(list(HUC12 = "150503010302", TOHUC = "150503010303",
                   comment =  "UNKNOWN fix (boader hack)")) %>%
    bind_rows(list(HUC12 = "111102010107", TOHUC = "111102010901",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111102010806", TOHUC = "111102010901",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111102010901", TOHUC = "111102010906",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111102010904", TOHUC = "111102010906",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111102010905", TOHUC = "111102010906",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "111102010906", TOHUC = "111102020401",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "100500010304", TOHUC = "100500010000",
                   comment =  "Milk Canada fix")) %>%
    bind_rows(list(HUC12 = "100500010000", TOHUC = "100500020000",
                   comment =  "Milk Canada fix")) %>%
    bind_rows(list(HUC12 = "100500020000", TOHUC = "100500020802",
                   comment =  "Milk Canada fix")) %>%
    bind_rows(list(HUC12 = "101303020905", TOHUC = "101303020908",
                   comment =  "Misdirected")) %>%
    bind_rows(list(HUC12 = "010801030205", TOHUC = "010801030206",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010801030204", TOHUC = "010801020403",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010801030203", TOHUC = "010801030205",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010801020403", TOHUC = "010801030205",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010801020402", TOHUC = "010801020403",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010801020401", TOHUC = "010801020403",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010801020302", TOHUC = "010801020403",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010801020203", TOHUC = "010801020403",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150701010202", TOHUC = "150701010203",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150701010109", TOHUC = "150701010204",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150701010203", TOHUC = "150701010204",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150701010204", TOHUC = "150701010207",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150702010804", TOHUC = "150702010805",
                   comment =  "misdirected")) %>%
    bind_rows(list(HUC12 = "090100060803", TOHUC = "090100060000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100060000", TOHUC = "090100070000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100070000", TOHUC = "090100080000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080000", TOHUC = "090100080703",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080703", TOHUC = "090100080705",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080804", TOHUC = "090100080807",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080807", TOHUC = "090100080810",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080810", TOHUC = "090100080905",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080905", TOHUC = "090100080908",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080908", TOHUC = "090100080909",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100080909", TOHUC = "090100081001",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "171100050302", TOHUC = "171100050303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010300010510", TOHUC = "010300010603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010300010603", TOHUC = "010300010610",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010300010610", TOHUC = "010300030101",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010300010609", TOHUC = "010300010610",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010300010608", TOHUC = "010300010610",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010300010604", TOHUC = "010300010610",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010300030101", TOHUC = "010300030106",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100090607", TOHUC = "CANADA",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100070404", TOHUC = "010100070405",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100070403", TOHUC = "010100070404",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "031002080204", TOHUC = "031002080205",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "031002080205", TOHUC = "031002080207",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080801020510", TOHUC = "080801020701",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051202081002", TOHUC = "051202081003",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051202081003", TOHUC = "051202081005",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051202080807", TOHUC = "051202081005",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051202081005", TOHUC = "051202081006",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "051202081001", TOHUC = "051202081003",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060101080205", TOHUC = "060101080604",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060101080603", TOHUC = "060101080604",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060101080601", TOHUC = "060101080604",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060101080306", TOHUC = "060101080604",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060101080604", TOHUC = "060101080606",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060102030107", TOHUC = "060102030303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060102030303", TOHUC = "060102030304",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060102030301", TOHUC = "060102030303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "060102030302", TOHUC = "060102030303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "180901030407", TOHUC = "180901030502",
                   comment =  "Owen Lake Partially Closed Basin")) %>%
    bind_rows(list(HUC12 = "031102060101", TOHUC = "031102060104",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "031102060103", TOHUC = "031102060104",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "031102060104", TOHUC = "031102060303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "031102060301", TOHUC = "031102060303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "031102060303", TOHUC = "031102060304",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170601070804", TOHUC = "170601070807",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170601070807", TOHUC = "170601070808",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170601070808", TOHUC = "170601100102",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170601070704", TOHUC = "170601100102",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170601070406", TOHUC = "170601070808",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170601070805", TOHUC = "170601070807",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170601070806", TOHUC = "170601070807",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100040203", TOHUC = "010100040603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100040603", TOHUC = "010100040606",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100040602", TOHUC = "010100040603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100040601", TOHUC = "010100040603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100040306", TOHUC = "010100040603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010100040104", TOHUC = "010100040603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302040705", TOHUC = "080302041104",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302041104", TOHUC = "080302020205",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302020205", TOHUC = "080302020600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302020600", TOHUC = "080302060104",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060106", TOHUC = "080302060300",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060300", TOHUC = "080302060407",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060407", TOHUC = "080302060703",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060706", TOHUC = "080302060802",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060802", TOHUC = "080302060805",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060805", TOHUC = "080302080103",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302080103", TOHUC = "080302080201",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060801", TOHUC = "080302060802",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060404", TOHUC = "080302060406",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060405", TOHUC = "080302060406",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060406", TOHUC = "080302060407",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060401", TOHUC = "080302060404",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060205", TOHUC = "080302060300",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302060900", TOHUC = "080302080103",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071700", TOHUC = "080302060900",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071900", TOHUC = "080302090600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302090600", TOHUC = "080302080201",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302090500", TOHUC = "080302090600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302090400", TOHUC = "080302090500",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302090300", TOHUC = "080302090400",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071600", TOHUC = "080302071700",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302090100", TOHUC = "080302090300",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302090200", TOHUC = "080302090300",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071400", TOHUC = "080302071600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071500", TOHUC = "080302071600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070900", TOHUC = "080302071000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071000", TOHUC = "080302071500",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070800", TOHUC = "080302071000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070500", TOHUC = "080302071000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071300", TOHUC = "080302071400",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071200", TOHUC = "080302071300",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071800", TOHUC = "080302090600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302071100", TOHUC = "080302071300",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070400", TOHUC = "080302070500",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070700", TOHUC = "080302070800",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070600", TOHUC = "080302070700",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070300", TOHUC = "080302070400",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070200", TOHUC = "080302070400",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302020400", TOHUC = "080302020600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302020500", TOHUC = "080302020600",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302070100", TOHUC = "080302070300",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302041000", TOHUC = "080302041104",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302040900", TOHUC = "080302041000",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "080302040800", TOHUC = "080302041104",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101302020610", TOHUC = "101302030101",
                   comment =  "Misdirected")) %>%
    bind_rows(list(HUC12 = "170200060101", TOHUC = "170200060108",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170200060108", TOHUC = "170200060211",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090203150208", TOHUC = "090203161003",
                   comment =  "CANADA connect")) %>%
    bind_rows(list(HUC12 = "160202030304", TOHUC = "160202030306",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "160202030306", TOHUC = "160202030405",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "160202030305", TOHUC = "160202030306",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "160202030302", TOHUC = "160202030306",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "160202030303", TOHUC = "160202030306",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150801010303", TOHUC = "150801010305",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150801010301", TOHUC = "150801010303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "150702010804", TOHUC = "150702010805",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "171100130204", TOHUC = "171100130205",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "171100130205", TOHUC = "171100130303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "171100130301", TOHUC = "171100130302",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "171100120104", TOHUC = "171100120105",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "171100120105", TOHUC = "171100120106",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "171100120106", TOHUC = "171100120400",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160103", TOHUC = "170102160105",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160104", TOHUC = "170102160105",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160105", TOHUC = "170102160107",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160107", TOHUC = "170102160111",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160101", TOHUC = "170102160107",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160102", TOHUC = "170102160107",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170102160106", TOHUC = "170102160107",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "101900030606", TOHUC = "101900030607",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170101060101", TOHUC = "170101010409",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170101010409", TOHUC = "170101010604",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170101010604", TOHUC = "170101010607",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "170101010607", TOHUC = "170101010708",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010400010203", TOHUC = "010400010602",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010400010406", TOHUC = "010400010602",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010400010502", TOHUC = "010400010602",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010400010601", TOHUC = "010400010602",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010400010602", TOHUC = "010400010603",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "010400010603", TOHUC = "010400010604",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100020301", TOHUC = "090100020303",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100060712", TOHUC = "090100060803",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100032203", TOHUC = "090100032207",
                   comment =  "UNKNOWN fix")) %>%
    bind_rows(list(HUC12 = "090100032207", TOHUC = "090100030000",
                   comment =  "UNKNOWN fix"))
  } else {
    stop("only 'latest' or 'nhdplusv2' version supported.")
  }
  return(fixes)
}


