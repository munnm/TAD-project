# This is based on depth first search

# An outer loop scans all nodes of graph and starts a search from every node

# Node neighbors are added to the cycle path

# Recursion ends if no more non-visited neighbors can be added

# A new cycle is found if the path is longer than 4 nodes and next neighbor is the start of the path

# To avoid duplicate cycles, the cycles are normalized by rotating the smallest node to the start



graph = [[1,2],[1,3],[1,4],[2,3],[3,4],[2,6],[4,6],[7,8],[8,9],[9,7]]
#d = 0.5  graph = [[1, 2], [1, 5], [1, 11], [1, 18], [1, 62], [1, 118], [3, 4], [3, 5], [6, 10], [6, 11], [7, 9], [7, 10], [8, 9], [12, 13], [12, 14], [12, 18], [15, 16], [15, 18], [17, 18], [19, 20], [19, 25], [19, 28], [19, 40], [19, 62], [21, 22], [21, 25], [23, 24], [23, 25], [26, 28], [27, 28], [29, 32], [29, 40], [30, 31], [30, 32], [33, 35], [33, 40], [34, 35], [36, 37], [36, 40], [38, 39], [38, 40], [41, 42], [41, 58], [41, 62], [43, 47], [43, 53], [43, 58], [44, 46], [44, 47], [45, 46], [48, 49], [48, 50], [48, 53], [51, 52], [51, 53], [54, 55], [54, 58], [56, 57], [56, 58], [59, 61], [59, 62], [60, 61], [63, 64], [63, 66], [63, 74], [63, 104], [63, 118], [65, 66], [67, 68], [67, 74], [69, 70], [69, 72], [69, 74], [71, 72], [73, 74], [75, 77], [75, 82], [75, 104], [76, 77], [78, 79], [78, 80], [78, 82], [81, 82], [83, 85], [83, 88], [83, 97], [83, 104], [84, 85], [86, 87], [86, 88], [89, 91], [89, 92], [89, 97], [90, 91], [93, 94], [93, 96], [93, 97], [95, 96], [98, 99], [98, 101], [98, 104], [100, 101], [102, 103], [102, 104], [105, 106], [105, 108], [105, 109], [105, 113], [105, 118], [107, 108], [110, 113], [111, 113], [112, 113], [114, 115], [114, 118], [116, 117], [116, 118], [119, 120], [119, 123], [119, 133], [119, 141], [119, 155], [119, 176], [121, 122], [121, 123], [124, 126], [124, 133], [125, 126], [127, 130], [127, 133], [128, 129], [128, 130], [131, 132], [131, 133], [134, 135], [134, 141], [136, 137], [136, 141], [138, 141], [139, 140], [139, 141], [142, 155], [143, 144], [143, 150], [143, 155], [145, 147], [145, 148], [145, 150], [146, 147], [149, 150], [151, 152], [151, 153], [151, 155], [154, 155], [156, 157], [156, 164], [156, 176], [158, 159], [158, 160], [158, 164], [161, 162], [161, 164], [163, 164], [165, 166], [165, 169], [165, 176], [167, 168], [167, 169], [170, 171], [170, 173], [170, 174], [170, 176], [172, 173], [175, 176], [177, 180], [177, 181], [177, 203], [177, 223], [177, 266], [178, 179], [178, 180], [182, 183], [182, 185], [182, 189], [182, 197], [182, 203], [184, 185], [186, 188], [186, 189], [187, 188], [190, 192], [190, 193], [190, 197], [191, 192], [194, 195], [194, 197], [196, 197], [198, 199], [198, 203], [200, 201], [200, 202], [200, 203], [204, 205], [204, 211], [204, 223], [206, 207], [206, 210], [206, 211], [208, 209], [208, 210], [212, 215], [212, 223], [213, 214], [213, 215], [216, 218], [216, 223], [217, 218], [219, 221], [219, 223], [220, 221], [222, 223], [224, 226], [224, 232], [224, 252], [224, 266], [225, 226], [227, 228], [227, 231], [227, 232], [229, 230], [229, 231], [233, 236], [233, 238], [233, 244], [233, 252], [234, 235], [234, 236], [237, 238], [239, 244], [240, 241], [240, 244], [242, 243], [242, 244], [245, 246], [245, 247], [245, 252], [248, 249], [248, 252], [250, 251], [250, 252], [253, 255], [253, 258], [253, 266], [254, 255], [256, 258], [257, 258], [259, 262], [259, 266], [260, 262], [261, 262], [263, 264], [263, 266], [265, 266], [267, 268], [267, 270], [267, 276], [267, 285], [267, 303], [267, 347], [267, 394], [269, 270], [271, 273], [271, 276], [272, 273], [274, 275], [274, 276], [277, 279], [277, 280], [277, 285], [278, 279], [281, 282], [281, 283], [281, 285], [284, 285], [286, 288], [286, 292], [286, 298], [286, 303], [287, 288], [289, 290], [289, 292], [291, 292], [293, 296], [293, 298], [294, 295], [294, 296], [297, 298], [299, 300], [299, 301], [299, 303], [302, 303], [304, 305], [304, 306], [304, 308], [304, 319], [304, 333], [304, 347], [307, 308], [309, 310], [309, 312], [309, 315], [309, 319], [311, 312], [313, 314], [313, 315], [316, 317], [316, 319], [318, 319], [320, 322], [320, 323], [320, 333], [321, 322], [324, 327], [324, 330], [324, 333], [325, 327], [326, 327], [328, 329], [328, 330], [331, 332], [331, 333], [334, 338], [334, 343], [334, 347], [335, 336], [335, 338], [337, 338], [339, 340], [339, 343], [341, 342], [341, 343], [344, 345], [344, 347], [346, 347], [348, 350], [348, 353], [348, 363], [348, 375], [348, 394], [349, 350], [351, 353], [352, 353], [354, 355], [354, 358], [354, 361], [354, 363], [356, 358], [357, 358], [359, 360], [359, 361], [362, 363], [364, 365], [364, 369], [364, 375], [366, 367], [366, 369], [368, 369], [370, 371], [370, 375], [372, 373], [372, 375], [374, 375], [376, 378], [376, 379], [376, 383], [376, 394], [377, 378], [380, 383], [381, 382], [381, 383], [384, 385], [384, 386], [384, 394], [387, 388], [387, 391], [387, 394], [389, 390], [389, 391], [392, 394], [393, 394]]
#d = 1    graph = [[1, 2], [1, 5], [1, 11], [1, 18], [1, 62], [1, 118], [3, 4], [3, 5], [6, 10], [6, 11], [7, 9], [7, 10], [8, 9], [12, 13], [12, 14], [12, 18], [15, 16], [15, 18], [17, 18], [19, 20], [19, 25], [19, 28], [19, 40], [19, 62], [21, 22], [21, 25], [23, 24], [23, 25], [26, 28], [27, 28], [29, 32], [29, 40], [30, 31], [30, 32], [33, 35], [33, 40], [34, 35], [36, 37], [36, 40], [38, 39], [38, 40], [41, 42], [41, 58], [41, 62], [43, 47], [43, 53], [43, 58], [44, 46], [44, 47], [45, 46], [48, 49], [48, 50], [48, 53], [51, 52], [51, 53], [54, 55], [54, 58], [56, 57], [56, 58], [59, 61], [59, 62], [60, 61], [63, 64], [63, 66], [63, 74], [63, 104], [63, 118], [65, 66], [67, 68], [67, 74], [69, 70], [69, 72], [69, 74], [71, 72], [73, 74], [75, 77], [75, 82], [75, 104], [76, 77], [78, 79], [78, 80], [78, 82], [81, 82], [83, 85], [83, 88], [83, 97], [83, 104], [84, 85], [86, 87], [86, 88], [89, 91], [89, 92], [89, 97], [90, 91], [93, 94], [93, 96], [93, 97], [95, 96], [98, 99], [98, 101], [98, 104], [100, 101], [102, 103], [102, 104], [105, 106], [105, 108], [105, 109], [105, 113], [105, 118], [107, 108], [110, 113], [111, 113], [112, 113], [114, 115], [114, 118], [116, 117], [116, 118], [119, 120], [119, 123], [119, 133], [119, 141], [119, 155], [119, 176], [119, 266], [121, 122], [121, 123], [124, 126], [124, 133], [125, 126], [127, 130], [127, 133], [128, 129], [128, 130], [131, 132], [131, 133], [134, 135], [134, 141], [136, 137], [136, 141], [138, 141], [139, 140], [139, 141], [142, 155], [143, 144], [143, 150], [143, 155], [145, 147], [145, 148], [145, 150], [146, 147], [149, 150], [151, 152], [151, 153], [151, 155], [154, 155], [156, 157], [156, 164], [156, 176], [158, 159], [158, 160], [158, 164], [161, 162], [161, 164], [163, 164], [165, 166], [165, 169], [165, 176], [167, 168], [167, 169], [170, 171], [170, 173], [170, 174], [170, 176], [172, 173], [175, 176], [177, 180], [177, 181], [177, 203], [177, 223], [177, 266], [178, 179], [178, 180], [182, 183], [182, 185], [182, 189], [182, 197], [182, 203], [184, 185], [186, 188], [186, 189], [187, 188], [190, 192], [190, 193], [190, 197], [191, 192], [194, 195], [194, 197], [196, 197], [198, 199], [198, 203], [200, 201], [200, 202], [200, 203], [204, 205], [204, 211], [204, 223], [206, 207], [206, 210], [206, 211], [208, 209], [208, 210], [212, 215], [212, 223], [213, 214], [213, 215], [216, 218], [216, 223], [217, 218], [219, 221], [219, 223], [220, 221], [222, 223], [224, 226], [224, 232], [224, 252], [224, 266], [225, 226], [227, 228], [227, 231], [227, 232], [229, 230], [229, 231], [233, 236], [233, 238], [233, 244], [233, 252], [234, 235], [234, 236], [237, 238], [239, 244], [240, 241], [240, 244], [242, 243], [242, 244], [245, 246], [245, 247], [245, 252], [248, 249], [248, 252], [250, 251], [250, 252], [253, 255], [253, 258], [253, 266], [254, 255], [256, 258], [257, 258], [259, 262], [259, 266], [260, 262], [261, 262], [263, 264], [263, 266], [265, 266], [267, 268], [267, 270], [267, 276], [267, 285], [267, 303], [267, 347], [267, 394], [269, 270], [271, 273], [271, 276], [272, 273], [274, 275], [274, 276], [277, 279], [277, 280], [277, 285], [278, 279], [281, 282], [281, 283], [281, 285], [284, 285], [286, 288], [286, 292], [286, 298], [286, 303], [287, 288], [289, 290], [289, 292], [291, 292], [293, 296], [293, 298], [294, 295], [294, 296], [297, 298], [299, 300], [299, 301], [299, 303], [302, 303], [304, 305], [304, 306], [304, 308], [304, 319], [304, 333], [304, 347], [307, 308], [309, 310], [309, 312], [309, 315], [309, 319], [311, 312], [313, 314], [313, 315], [316, 317], [316, 319], [318, 319], [320, 322], [320, 323], [320, 333], [321, 322], [324, 327], [324, 330], [324, 333], [325, 327], [326, 327], [328, 329], [328, 330], [331, 332], [331, 333], [334, 338], [334, 343], [334, 347], [335, 336], [335, 338], [337, 338], [339, 340], [339, 343], [341, 342], [341, 343], [344, 345], [344, 347], [346, 347], [348, 350], [348, 353], [348, 363], [348, 375], [348, 394], [349, 350], [351, 353], [352, 353], [354, 355], [354, 358], [354, 361], [354, 363], [356, 358], [357, 358], [359, 360], [359, 361], [362, 363], [364, 365], [364, 369], [364, 375], [366, 367], [366, 369], [368, 369], [370, 371], [370, 375], [372, 373], [372, 375], [374, 375], [376, 378], [376, 379], [376, 383], [376, 394], [377, 378], [380, 383], [381, 382], [381, 383], [384, 385], [384, 386], [384, 394], [387, 388], [387, 391], [387, 394], [389, 390], [389, 391], [392, 394], [393, 394]]
#d = 2.5  graph = [[1, 2], [1, 5], [1, 11], [1, 18], [1, 62], [1, 118], [3, 4], [3, 5], [6, 10], [6, 11], [7, 9], [7, 10], [8, 9], [12, 13], [12, 14], [12, 18], [15, 16], [15, 18], [17, 18], [19, 20], [19, 25], [19, 28], [19, 40], [19, 62], [21, 22], [21, 25], [23, 24], [23, 25], [26, 28], [27, 28], [29, 32], [29, 40], [30, 31], [30, 32], [33, 35], [33, 40], [34, 35], [36, 37], [36, 40], [38, 39], [38, 40], [41, 42], [41, 58], [41, 62], [43, 47], [43, 53], [43, 58], [44, 46], [44, 47], [45, 46], [48, 49], [48, 50], [48, 53], [51, 52], [51, 53], [54, 55], [54, 58], [56, 57], [56, 58], [59, 61], [59, 62], [60, 61], [63, 64], [63, 66], [63, 74], [63, 104], [63, 118], [65, 66], [67, 68], [67, 74], [69, 70], [69, 72], [69, 74], [71, 72], [73, 74], [75, 77], [75, 82], [75, 104], [76, 77], [78, 79], [78, 80], [78, 82], [81, 82], [83, 85], [83, 88], [83, 97], [83, 104], [84, 85], [86, 87], [86, 88], [89, 91], [89, 92], [89, 97], [90, 91], [93, 94], [93, 96], [93, 97], [95, 96], [98, 99], [98, 101], [98, 104], [100, 101], [102, 103], [102, 104], [105, 106], [105, 108], [105, 109], [105, 113], [105, 118], [107, 108], [110, 113], [111, 113], [112, 113], [114, 115], [114, 118], [116, 117], [116, 118], [119, 120], [119, 123], [119, 133], [119, 141], [119, 155], [119, 176], [119, 266], [119, 394], [121, 122], [121, 123], [124, 126], [124, 133], [125, 126], [127, 130], [127, 133], [128, 129], [128, 130], [131, 132], [131, 133], [134, 135], [134, 141], [136, 137], [136, 141], [138, 141], [139, 140], [139, 141], [142, 155], [143, 144], [143, 150], [143, 155], [145, 147], [145, 148], [145, 150], [146, 147], [149, 150], [151, 152], [151, 153], [151, 155], [154, 155], [156, 157], [156, 164], [156, 176], [158, 159], [158, 160], [158, 164], [161, 162], [161, 164], [163, 164], [165, 166], [165, 169], [165, 176], [167, 168], [167, 169], [170, 171], [170, 173], [170, 174], [170, 176], [172, 173], [175, 176], [177, 180], [177, 181], [177, 203], [177, 223], [177, 266], [178, 179], [178, 180], [182, 183], [182, 185], [182, 189], [182, 197], [182, 203], [184, 185], [186, 188], [186, 189], [187, 188], [190, 192], [190, 193], [190, 197], [191, 192], [194, 195], [194, 197], [196, 197], [198, 199], [198, 203], [200, 201], [200, 202], [200, 203], [204, 205], [204, 211], [204, 223], [206, 207], [206, 210], [206, 211], [208, 209], [208, 210], [212, 215], [212, 223], [213, 214], [213, 215], [216, 218], [216, 223], [217, 218], [219, 221], [219, 223], [220, 221], [222, 223], [224, 226], [224, 232], [224, 252], [224, 266], [225, 226], [227, 228], [227, 231], [227, 232], [229, 230], [229, 231], [233, 236], [233, 238], [233, 244], [233, 252], [234, 235], [234, 236], [237, 238], [239, 244], [240, 241], [240, 244], [242, 243], [242, 244], [245, 246], [245, 247], [245, 252], [248, 249], [248, 252], [250, 251], [250, 252], [253, 255], [253, 258], [253, 266], [254, 255], [256, 258], [257, 258], [259, 262], [259, 266], [260, 262], [261, 262], [263, 264], [263, 266], [265, 266], [267, 268], [267, 270], [267, 276], [267, 285], [267, 303], [267, 347], [267, 394], [269, 270], [271, 273], [271, 276], [272, 273], [274, 275], [274, 276], [277, 279], [277, 280], [277, 285], [278, 279], [281, 282], [281, 283], [281, 285], [284, 285], [286, 288], [286, 292], [286, 298], [286, 303], [287, 288], [289, 290], [289, 292], [291, 292], [293, 296], [293, 298], [294, 295], [294, 296], [297, 298], [299, 300], [299, 301], [299, 303], [302, 303], [304, 305], [304, 306], [304, 308], [304, 319], [304, 333], [304, 347], [307, 308], [309, 310], [309, 312], [309, 315], [309, 319], [311, 312], [313, 314], [313, 315], [316, 317], [316, 319], [318, 319], [320, 322], [320, 323], [320, 333], [321, 322], [324, 327], [324, 330], [324, 333], [325, 327], [326, 327], [328, 329], [328, 330], [331, 332], [331, 333], [334, 338], [334, 343], [334, 347], [335, 336], [335, 338], [337, 338], [339, 340], [339, 343], [341, 342], [341, 343], [344, 345], [344, 347], [346, 347], [348, 350], [348, 353], [348, 363], [348, 375], [348, 394], [349, 350], [351, 353], [352, 353], [354, 355], [354, 358], [354, 361], [354, 363], [356, 358], [357, 358], [359, 360], [359, 361], [362, 363], [364, 365], [364, 369], [364, 375], [366, 367], [366, 369], [368, 369], [370, 371], [370, 375], [372, 373], [372, 375], [374, 375], [376, 378], [376, 379], [376, 383], [376, 394], [377, 378], [380, 383], [381, 382], [381, 383], [384, 385], [384, 386], [384, 394], [387, 388], [387, 391], [387, 394], [389, 390], [389, 391], [392, 394], [393, 394]]
#d = 3    graph = [[1, 2], [1, 5], [1, 11], [1, 18], [1, 62], [1, 118], [1, 394], [3, 4], [3, 5], [6, 10], [6, 11], [7, 9], [7, 10], [8, 9], [12, 13], [12, 14], [12, 18], [15, 16], [15, 18], [17, 18], [19, 20], [19, 25], [19, 28], [19, 40], [19, 62], [21, 22], [21, 25], [23, 24], [23, 25], [26, 28], [27, 28], [29, 32], [29, 40], [30, 31], [30, 32], [33, 35], [33, 40], [34, 35], [36, 37], [36, 40], [38, 39], [38, 40], [41, 42], [41, 58], [41, 62], [43, 47], [43, 53], [43, 58], [44, 46], [44, 47], [45, 46], [48, 49], [48, 50], [48, 53], [51, 52], [51, 53], [54, 55], [54, 58], [56, 57], [56, 58], [59, 61], [59, 62], [60, 61], [63, 64], [63, 66], [63, 74], [63, 104], [63, 118], [65, 66], [67, 68], [67, 74], [69, 70], [69, 72], [69, 74], [71, 72], [73, 74], [75, 77], [75, 82], [75, 104], [76, 77], [78, 79], [78, 80], [78, 82], [81, 82], [83, 85], [83, 88], [83, 97], [83, 104], [84, 85], [86, 87], [86, 88], [89, 91], [89, 92], [89, 97], [90, 91], [93, 94], [93, 96], [93, 97], [95, 96], [98, 99], [98, 101], [98, 104], [100, 101], [102, 103], [102, 104], [105, 106], [105, 108], [105, 109], [105, 113], [105, 118], [107, 108], [110, 113], [111, 113], [112, 113], [114, 115], [114, 118], [116, 117], [116, 118], [119, 120], [119, 123], [119, 133], [119, 141], [119, 155], [119, 176], [119, 266], [119, 394], [121, 122], [121, 123], [124, 126], [124, 133], [125, 126], [127, 130], [127, 133], [128, 129], [128, 130], [131, 132], [131, 133], [134, 135], [134, 141], [136, 137], [136, 141], [138, 141], [139, 140], [139, 141], [142, 155], [143, 144], [143, 150], [143, 155], [145, 147], [145, 148], [145, 150], [146, 147], [149, 150], [151, 152], [151, 153], [151, 155], [154, 155], [156, 157], [156, 164], [156, 176], [158, 159], [158, 160], [158, 164], [161, 162], [161, 164], [163, 164], [165, 166], [165, 169], [165, 176], [167, 168], [167, 169], [170, 171], [170, 173], [170, 174], [170, 176], [172, 173], [175, 176], [177, 180], [177, 181], [177, 203], [177, 223], [177, 266], [178, 179], [178, 180], [182, 183], [182, 185], [182, 189], [182, 197], [182, 203], [184, 185], [186, 188], [186, 189], [187, 188], [190, 192], [190, 193], [190, 197], [191, 192], [194, 195], [194, 197], [196, 197], [198, 199], [198, 203], [200, 201], [200, 202], [200, 203], [204, 205], [204, 211], [204, 223], [206, 207], [206, 210], [206, 211], [208, 209], [208, 210], [212, 215], [212, 223], [213, 214], [213, 215], [216, 218], [216, 223], [217, 218], [219, 221], [219, 223], [220, 221], [222, 223], [224, 226], [224, 232], [224, 252], [224, 266], [225, 226], [227, 228], [227, 231], [227, 232], [229, 230], [229, 231], [233, 236], [233, 238], [233, 244], [233, 252], [234, 235], [234, 236], [237, 238], [239, 244], [240, 241], [240, 244], [242, 243], [242, 244], [245, 246], [245, 247], [245, 252], [248, 249], [248, 252], [250, 251], [250, 252], [253, 255], [253, 258], [253, 266], [254, 255], [256, 258], [257, 258], [259, 262], [259, 266], [260, 262], [261, 262], [263, 264], [263, 266], [265, 266], [267, 268], [267, 270], [267, 276], [267, 285], [267, 303], [267, 347], [267, 394], [269, 270], [271, 273], [271, 276], [272, 273], [274, 275], [274, 276], [277, 279], [277, 280], [277, 285], [278, 279], [281, 282], [281, 283], [281, 285], [284, 285], [286, 288], [286, 292], [286, 298], [286, 303], [287, 288], [289, 290], [289, 292], [291, 292], [293, 296], [293, 298], [294, 295], [294, 296], [297, 298], [299, 300], [299, 301], [299, 303], [302, 303], [304, 305], [304, 306], [304, 308], [304, 319], [304, 333], [304, 347], [307, 308], [309, 310], [309, 312], [309, 315], [309, 319], [311, 312], [313, 314], [313, 315], [316, 317], [316, 319], [318, 319], [320, 322], [320, 323], [320, 333], [321, 322], [324, 327], [324, 330], [324, 333], [325, 327], [326, 327], [328, 329], [328, 330], [331, 332], [331, 333], [334, 338], [334, 343], [334, 347], [335, 336], [335, 338], [337, 338], [339, 340], [339, 343], [341, 342], [341, 343], [344, 345], [344, 347], [346, 347], [348, 350], [348, 353], [348, 363], [348, 375], [348, 394], [349, 350], [351, 353], [352, 353], [354, 355], [354, 358], [354, 361], [354, 363], [356, 358], [357, 358], [359, 360], [359, 361], [362, 363], [364, 365], [364, 369], [364, 375], [366, 367], [366, 369], [368, 369], [370, 371], [370, 375], [372, 373], [372, 375], [374, 375], [376, 378], [376, 379], [376, 383], [376, 394], [377, 378], [380, 383], [381, 382], [381, 383], [384, 385], [384, 386], [384, 394], [387, 388], [387, 391], [387, 394], [389, 390], [389, 391], [392, 394], [393, 394]]
#d= 0.5, Neuron
#graph = [[1, 2], [1, 4], [1, 12], [1, 25], [1, 55], [3, 4], [5, 7], [5, 8], [5, 12], [6, 7], [9, 10], [9, 12], [11, 12], [13, 14], [13, 16], [13, 25], [15, 16], [17, 18], [17, 25], [19, 24], [19, 25], [20, 23], [20, 24], [21, 23], [22, 23], [26, 27], [26, 29], [26, 36], [26, 55], [28, 29], [30, 32], [30, 35], [30, 36], [31, 32], [33, 35], [34, 35], [37, 40], [37, 43], [37, 50], [37, 55], [38, 39], [38, 40], [41, 42], [41, 43], [44, 46], [44, 47], [44, 49], [44, 50], [45, 46], [48, 49], [51, 52], [51, 55], [53, 54], [53, 55], [56, 60], [56, 62], [56, 76], [56, 114], [57, 60], [58, 60], [59, 60], [61, 62], [63, 64], [63, 68], [63, 69], [63, 72], [63, 76], [65, 68], [66, 68], [67, 68], [70, 72], [71, 72], [73, 74], [73, 76], [75, 76], [77, 79], [77, 103], [77, 114], [78, 79], [80, 81], [80, 84], [80, 86], [80, 101], [80, 103], [82, 83], [82, 84], [85, 86], [87, 88], [87, 93], [87, 101], [89, 91], [89, 92], [89, 93], [90, 91], [94, 95], [94, 96], [94, 100], [94, 101], [97, 98], [97, 100], [99, 100], [102, 103], [104, 105], [104, 109], [104, 114], [106, 107], [106, 108], [106, 109], [110, 111], [110, 114], [112, 113], [112, 114], [115, 117], [115, 118], [115, 121], [115, 127], [115, 144], [116, 117], [119, 120], [119, 121], [122, 123], [122, 127], [124, 126], [124, 127], [125, 126], [128, 129], [128, 132], [128, 134], [128, 144], [130, 132], [131, 132], [133, 134], [135, 136], [135, 137], [135, 144], [138, 141], [138, 144], [139, 141], [140, 141], [142, 144], [143, 144], [145, 147], [145, 157], [146, 147], [148, 149], [148, 150], [148, 152], [148, 153], [148, 154], [148, 157], [151, 152], [155, 156], [155, 157], [158, 159], [158, 160], [158, 165], [158, 177], [158, 265], [161, 162], [161, 164], [161, 165], [163, 164], [166, 167], [166, 169], [166, 172], [166, 177], [168, 169], [170, 172], [171, 172], [173, 174], [173, 177], [175, 176], [175, 177], [178, 180], [178, 186], [178, 209], [178, 240], [178, 265], [179, 180], [181, 184], [181, 186], [182, 183], [182, 184], [185, 186], [187, 188], [187, 193], [187, 198], [187, 203], [187, 209], [189, 190], [189, 193], [191, 193], [192, 193], [194, 195], [194, 196], [194, 198], [197, 198], [199, 200], [199, 202], [199, 203], [201, 202], [204, 205], [204, 206], [204, 209], [207, 208], [207, 209], [210, 212], [210, 224], [210, 225], [210, 240], [211, 212], [213, 215], [213, 221], [213, 224], [214, 215], [216, 219], [216, 221], [217, 219], [218, 219], [220, 221], [222, 223], [222, 224], [226, 231], [226, 236], [226, 240], [227, 229], [227, 231], [228, 229], [230, 231], [232, 235], [232, 236], [233, 235], [234, 235], [237, 238], [237, 240], [239, 240], [241, 243], [241, 246], [241, 265], [242, 243], [244, 245], [244, 246], [247, 248], [247, 253], [247, 256], [247, 265], [249, 252], [249, 253], [250, 252], [251, 252], [254, 255], [254, 256], [257, 258], [257, 261], [257, 265], [259, 261], [260, 261], [262, 263], [262, 265], [264, 265], [266, 270], [266, 285], [266, 334], [267, 268], [267, 270], [269, 270], [271, 272], [271, 274], [271, 275], [271, 280], [271, 285], [273, 274], [276, 277], [276, 280], [278, 279], [278, 280], [281, 282], [281, 284], [281, 285], [283, 284], [286, 287], [286, 288], [286, 290], [286, 298], [286, 319], [286, 334], [289, 290], [291, 292], [291, 293], [291, 298], [294, 295], [294, 298], [296, 297], [296, 298], [299, 301], [299, 319], [300, 301], [302, 303], [302, 308], [302, 319], [304, 305], [304, 306], [304, 308], [307, 308], [309, 310], [309, 312], [309, 318], [309, 319], [311, 312], [313, 314], [313, 318], [315, 316], [315, 318], [317, 318], [320, 323], [320, 327], [320, 334], [321, 323], [322, 323], [324, 327], [325, 326], [325, 327], [328, 329], [328, 330], [328, 334], [331, 332], [331, 333], [331, 334], [335, 336], [335, 340], [335, 343], [335, 360], [335, 371], [337, 338], [337, 339], [337, 340], [341, 342], [341, 343], [344, 347], [344, 354], [344, 360], [345, 346], [345, 347], [348, 350], [348, 354], [349, 350], [351, 354], [352, 353], [352, 354], [355, 360], [356, 357], [356, 360], [358, 359], [358, 360], [361, 362], [361, 364], [361, 367], [361, 371], [363, 364], [365, 366], [365, 367], [368, 369], [368, 371], [370, 371], [372, 374], [372, 379], [372, 386], [372, 402], [372, 415], [373, 374], [375, 376], [375, 379], [377, 378], [377, 379], [380, 382], [380, 386], [381, 382], [383, 384], [383, 386], [385, 386], [387, 388], [387, 389], [387, 394], [387, 399], [387, 402], [390, 391], [390, 392], [390, 394], [393, 394], [395, 396], [395, 399], [397, 398], [397, 399], [400, 401], [400, 402], [403, 405], [403, 408], [403, 415], [404, 405], [406, 407], [406, 408], [409, 410], [409, 412], [409, 415], [411, 412], [413, 415], [414, 415]]
#graph = [[1, 2], [1, 4], [1, 12], [1, 25], [3, 4], [5, 7], [5, 8], [5, 12], [6, 7], [9, 10], [9, 12], [11, 12], [13, 14], [13, 16], [13, 25], [15, 16], [17, 18], [17, 25], [19, 24], [19, 25], [20, 23], [20, 24], [21, 23], [22, 23], [26, 27], [26, 29], [26, 36], [26, 55], [28, 29], [30, 32], [30, 35], [30, 36], [31, 32], [33, 35], [34, 35], [37, 40], [37, 43], [37, 50], [37, 55], [38, 39], [38, 40], [41, 42], [41, 43], [44, 46], [44, 47], [44, 49], [44, 50], [45, 46], [48, 49], [51, 52], [51, 55], [53, 54], [53, 55], [56, 60], [56, 62], [57, 60], [58, 60], [59, 60], [61, 62], [63, 64], [63, 68], [63, 69], [63, 72], [63, 76], [65, 68], [66, 68], [67, 68], [70, 72], [71, 72], [73, 74], [73, 76], [75, 76], [77, 79], [77, 103], [77, 114], [78, 79], [80, 81], [80, 84], [80, 86], [80, 101], [80, 103], [82, 83], [82, 84], [85, 86], [87, 88], [87, 93], [87, 101], [89, 91], [89, 92], [89, 93], [90, 91], [94, 95], [94, 96], [94, 100], [94, 101], [97, 98], [97, 100], [99, 100], [102, 103], [104, 105], [104, 109], [104, 114], [106, 107], [106, 108], [106, 109], [110, 111], [110, 114], [112, 113], [112, 114], [115, 117], [115, 118], [115, 121], [115, 127], [116, 117], [119, 120], [119, 121], [122, 123], [122, 127], [124, 126], [124, 127], [125, 126], [128, 129], [128, 132], [128, 134], [128, 144], [130, 132], [131, 132], [133, 134], [135, 136], [135, 137], [135, 144], [138, 141], [138, 144], [139, 141], [140, 141], [142, 144], [143, 144], [145, 147], [145, 157], [146, 147], [148, 149], [148, 150], [148, 152], [148, 153], [148, 154], [148, 157], [151, 152], [155, 156], [155, 157], [158, 159], [158, 160], [158, 165], [158, 177], [161, 162], [161, 164], [161, 165], [163, 164], [166, 167], [166, 169], [166, 172], [166, 177], [168, 169], [170, 172], [171, 172], [173, 174], [173, 177], [175, 176], [175, 177], [178, 180], [178, 186], [178, 209], [179, 180], [181, 184], [181, 186], [182, 183], [182, 184], [185, 186], [187, 188], [187, 193], [187, 198], [187, 203], [187, 209], [189, 190], [189, 193], [191, 193], [192, 193], [194, 195], [194, 196], [194, 198], [197, 198], [199, 200], [199, 202], [199, 203], [201, 202], [204, 205], [204, 206], [204, 209], [207, 208], [207, 209], [210, 212], [210, 224], [210, 225], [210, 240], [211, 212], [213, 215], [213, 221], [213, 224], [214, 215], [216, 219], [216, 221], [217, 219], [218, 219], [220, 221], [222, 223], [222, 224], [226, 231], [226, 236], [226, 240], [227, 229], [227, 231], [228, 229], [230, 231], [232, 235], [232, 236], [233, 235], [234, 235], [237, 238], [237, 240], [239, 240], [241, 243], [241, 246], [242, 243], [244, 245], [244, 246], [247, 248], [247, 253], [247, 256], [247, 265], [249, 252], [249, 253], [250, 252], [251, 252], [254, 255], [254, 256], [257, 258], [257, 261], [257, 265], [259, 261], [260, 261], [262, 263], [262, 265], [264, 265], [266, 270], [266, 285], [267, 268], [267, 270], [269, 270], [271, 272], [271, 274], [271, 275], [271, 280], [271, 285], [273, 274], [276, 277], [276, 280], [278, 279], [278, 280], [281, 282], [281, 284], [281, 285], [283, 284], [286, 287], [286, 288], [286, 290], [286, 298], [286, 319], [289, 290], [291, 292], [291, 293], [291, 298], [294, 295], [294, 298], [296, 297], [296, 298], [299, 301], [299, 319], [300, 301], [302, 303], [302, 308], [302, 319], [304, 305], [304, 306], [304, 308], [307, 308], [309, 310], [309, 312], [309, 318], [309, 319], [311, 312], [313, 314], [313, 318], [315, 316], [315, 318], [317, 318], [320, 323], [320, 327], [320, 334], [321, 323], [322, 323], [324, 327], [325, 326], [325, 327], [328, 329], [328, 330], [328, 334], [331, 332], [331, 333], [331, 334], [335, 336], [335, 340], [335, 343], [335, 360], [337, 338], [337, 339], [337, 340], [341, 342], [341, 343], [344, 347], [344, 354], [344, 360], [345, 346], [345, 347], [348, 350], [348, 354], [349, 350], [351, 354], [352, 353], [352, 354], [355, 360], [356, 357], [356, 360], [358, 359], [358, 360], [361, 362], [361, 364], [361, 367], [361, 371], [363, 364], [365, 366], [365, 367], [368, 369], [368, 371], [370, 371], [372, 374], [372, 379], [372, 386], [373, 374], [375, 376], [375, 379], [377, 378], [377, 379], [380, 382], [380, 386], [381, 382], [383, 384], [383, 386], [385, 386], [387, 388], [387, 389], [387, 394], [387, 399], [387, 402], [390, 391], [390, 392], [390, 394], [393, 394], [395, 396], [395, 399], [397, 398], [397, 399], [400, 401], [400, 402], [403, 405], [403, 408], [403, 415], [404, 405], [406, 407], [406, 408], [409, 410], [409, 412], [409, 415], [411, 412], [413, 415], [414, 415]]
# can try testing this out on differing values for d
cycles = []


def main():
    global graph
    global cycles
    # first sort graph so that node1 < node2 for all edges in graph
    for edge in graph:
        edge = edge.sort()
    # order edges in graph by ascending order on node1
    graph = sorted(graph)    
    for edge in graph:
      for vertex in edge:
        findNewCycles([vertex])
    print(cycles)

def findNewCycles(path):
    start_node = path[0]
    next_node = None
    sub = []

    for edge in graph:
        node1, node2 = edge
        # stop search once you're too far along in the graph
        if node1 <= path[0]:
            if start_node in edge:
                    if node1 == start_node:
                        next_node = node2
                    else:
                        next_node = node1
            if not visited(next_node, path):
		# neighbor node not on path yet
                    sub = [next_node]
                    sub.extend(path)
		# explore extended path
                    findNewCycles(sub);
            elif len(path) > 3  and next_node == path[-1]:
		# cycle found
                    p = rotate_to_smallest(path);
                    inv = invert(p)
                    if isNew(p) and isNew(inv):
                        cycles.append(p)


# next two function is for formating purpose only
def invert(path):
    return rotate_to_smallest(path[::-1])



# rotate cycle path such that it begins with the smallest node
def rotate_to_smallest(path):
    n = path.index(min(path))
    return path[n:]+path[:n]



# check if it is a unique path
def isNew(path):
    return (not path in cycles)



# check if a vertex has been visited
def visited(vertex, path):
    return (vertex in path)

main()