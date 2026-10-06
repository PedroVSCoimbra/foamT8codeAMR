# Integração t8code ↔ foam-extend-5.0

## 1. Objetivo

Integrar o **t8code** como backend de malha adaptativa no
foam-extend-5.0, substituindo o motor `hexRef8` por um sistema
multi-elemento, paralelo e com load balancing.

## 2. Arquitetura

### Hierarquia de classes

```text
topoChangerFvMesh (foam-extend — base para load balancing)
  └── dynamicT8codeFvMesh (nova — t8code como backend completo)
```

Herdar de `topoChangerFvMesh` reutiliza a infraestrutura de
redistribuição paralela já existente no foam-extend.

### Princípio de design

O t8code segue o *Hollywood principle*: a aplicação registra
callbacks e o t8code os invoca durante Adapt, Interpolate,
Partition, Ghost e Iterate. A integração usa esse mecanismo
nativo em vez de reimplementar operações manualmente.

### Compatibilidade com OpenFOAM

O `fvMesh::updateMesh(mapPolyMesh)` faz quatro operações
que não podem ser bypassed (`fvMesh.C:627-655`):

```cpp
polyMesh::updateMesh(mpm);         // atualiza topologia
surfaceInterpolation::clearOut();  // invalida geometria
mapFields(mpm);                    // mapeia campos
mapOldVolumes(mpm);                // volumes para ddt
```

A estratégia é **híbrida**: t8code gerencia a forest e os
dados via callbacks; a conversão forest → `directTopoChange`
produz o `mapPolyMesh` para compatibilidade com
`fvMesh::updateMesh()`.

```text
t8code (forest + callbacks)
  → directTopoChange (compatibilidade)
    → mapPolyMesh
      → fvMesh::updateMesh
```

## 3. Fluxo por passo de tempo

```text
interDyMFoam.C
  └── mesh.update()  ← dynamicT8codeFvMesh::update()
        │
        ├── 1. Copiar campos OpenFOAM → arrays t8code
        │     (alpha1, U, pd como arrays SFC-ordered)
        │     t8_forest_set_user_data(forest, &app_data)
        │
        ├── 2. Adapt (t8code chama callback)
        │     t8_forest_set_adapt() + t8_forest_commit()
        │     Callback decide: 1=refina, -1=coarsen, 0=mantém
        │
        ├── 3. Balance (2:1)
        │     t8_forest_set_balance() + t8_forest_commit()
        │
        ├── 4. Interpolate (t8code chama callback)
        │     Mapeia campos entre malha velha e nova
        │     Refino: valor_pai + grad·offset
        │     Coarsen: média volumétrica
        │
        ├── 5. Partition + Partition Data
        │     t8_forest_set_partition() + t8_forest_commit()
        │     Arrays de dados redistribuídos automaticamente
        │
        ├── 6. Ghost + Ghost exchange
        │     t8_forest_set_ghost() + t8_forest_commit()
        │     Dados sincronizados entre processos
        │
        ├── 7. Converter forest → directTopoChange
        │     t8_forest_iterate() constrói points, faces,
        │     owner, neighbour via callback
        │     → polyAddPoint/Face/Cell no directTopoChange
        │
        ├── 8. fvMesh::updateMesh(map)  (OpenFOAM)
        │     mapFields, mapOldVolumes, clearOut
        │
        └── 9. Solver continua (interDyMFoam)
              correctPhi.H, gh, meshCourantNo, PIMPLE
```

## 4. Callbacks da integração

### Adapt

```c
int adapt_cb(forest, forest_from, which_tree, tree_class,
             lelement_id, scheme, is_family, num_elem, elements[])
{
    AppData *d = t8_forest_get_user_data(forest);
    double a1 = d->alpha1[lelement_id];

    if (is_family && should_coarsen(a1, num_elem))
        return -1;
    if (should_refine(a1))
        return 1;
    return 0;
}
```

### Interpolate (replace)

```c
void replace_cb(forest_old, forest_new, which_tree,
                num_outgoing, first_outgoing,
                num_incoming, first_incoming, user_data)
{
    // Refino (1→N): valor_pai + grad·(centro_filho - centro_pai)
    // Coarsen (N→1): soma(vol_i * valor_i) / soma(vol_i)
    // t8_forest_element_volume() dá volumes exatos
}
```

### Iterate (construção da malha OpenFOAM)

```c
void iterate_cb(forest, user_data, which_tree, tree_class,
                lelement_id, is_local, is_ghost, elem)
{
    // t8_forest_element_coordinate() → pontos
    // t8_element_get_face_corner()   → faces
    // t8_forest_element_face_neighbor() → owner/neighbour
    //   -1 = fronteira, -2 = processor patch
}
```

### User data

```c
typedef struct {
    double *alpha1;  // array[n_elements], SFC-ordered
    double *U;       // array[3*n_elements]
    double *pd;      // array[n_elements]
} AppData;

AppData data;
data.alpha1 = malloc(n * sizeof(double));
t8_forest_set_user_data(forest, &data);
```

## 5. API t8code (verificada no código-fonte)

| Operação | Função | Header |
| --- | --- | --- |
| Init forest | `t8_forest_init(&forest)` | `t8_forest_general.h:139` |
| Coarse mesh | `t8_forest_set_cmesh(f, cmesh, comm)` | `t8_forest_general.h:198` |
| Adapt | `t8_forest_set_adapt(f, from, cb, recursive)` | `t8_forest_general.h:266` |
| Balance 2:1 | `t8_forest_set_balance(f, from, repart)` | `t8_forest_general.h:368` |
| Partition | `t8_forest_set_partition(f, from, flag)` | `t8_forest_general.h:332` |
| Ghost | `t8_forest_set_ghost(f, type)` | `t8_forest_general.h` |
| Commit | `t8_forest_commit(f)` | `t8_forest_general.h:419` |
| Nº elem locais | `t8_forest_get_local_num_leaf_elements(f)` | `t8_forest_general.h:441` |
| Obter elemento | `t8_forest_get_leaf_element(f, id, &tree)` | `t8_forest_general.h:823` |
| Coordenada vértice | `t8_forest_element_coordinate(f, t, e, c, coords)` | `t8_forest_geometrical.h:54` |
| Vizinho de face | `t8_forest_element_face_neighbor(...)` | `t8_forest_general.h:923` |
| Volume elemento | `t8_forest_element_volume(f, t, e)` | `t8_forest_geometrical.h:144` |
| User data | `t8_forest_set/get_user_data(f, ptr)` | `t8_forest_general.h:276,285` |
| Iterate | `t8_forest_iterate(f, data, cb)` | `t8_forest_iterate.h` |
| Replace | `t8_forest_iterate_replace(f, cb)` | `t8_forest_iterate.h:212` |

## 6. Mapeamento de campos

### Refino (1→N): interpolação linear reconstruída

```text
valor_filho = valor_pai + grad(alpha1) · (centro_filho - centro_pai)
valor_filho = clamp(valor_filho, 0, 1)
```

O gradiente é estimado via ghost layer do t8code. Isso reduz o
alargamento da interface VOF.

### Coarsening (N→1): média volumétrica

```text
valor_pai = soma(vol_i · valor_i) / soma(vol_i)
```

`t8_forest_element_volume()` fornece volumes exatos.
Conservação `∫ α₁ dV` garantida geometricamente.

## 7. Paralelismo

O t8code gerencia todo o paralelismo MPI nativamente:

- **Partition:** reparticiona por SFC (escalável para 1M
  processos)
- **Partition Data:** redistribui arrays de dados junto
  com a malha
- **Ghost:** computa ghost layer entre processos
- **Ghost exchange:** sincroniza dados entre vizinhos MPI

O OpenFOAM apenas converte a forest local para
`directTopoChange` e reconstrói processor patches nas
interfaces (quando `t8_forest_element_face_neighbor` retorna
-2).

## 8. Conversão de malha

### Início: polyMesh → Forest

Apenas hexaedros puros (blockMesh) na Fase 1. Cada célula
hexa vira um elemento `T8_ECLASS_HEX` na coarse mesh.

### Cada passo: Forest → directTopoChange

Via `t8_forest_iterate()`, o callback extrai para cada
elemento:

1. Vértices (`t8_forest_element_coordinate`)
2. Faces (`t8_element_get_face_corner`)
3. Vizinhos (`t8_forest_element_face_neighbor`)
4. Patches de contorno (neighbor = -1 → wall/symmetry,
   neighbor = -2 → processor)

## 9. Dicionário de configuração

```c++
dynamicFvMesh   dynamicT8codeFvMesh;

libs            ("libdynamicT8codeFvMesh.so");

dynamicT8codeFvMeshCoeffs
{
    refineInterval      1;
    field               alpha1;
    lowerRefineLevel    0.001;
    upperRefineLevel    0.999;
    unrefineLevel       10;
    nBufferLayers       1;
    maxRefinement       2;
    maxCells            200000;
    correctFluxes       ((phi U));
    dumpLevel           true;

    elementType         hexahedron;
    balanceType         2:1;
    partitionEnable     true;
    imbalanceThreshold  0.8;
    ghostType           face;
}
```

## 10. Arquivos

### Novos

```text
src/dynamicMesh/topoChangerFvMesh/dynamicT8codeFvMesh/
├── dynamicT8codeFvMesh.H
├── dynamicT8codeFvMesh.C
├── T8codeCallbacks.H
├── T8codeCallbacks.C
├── T8codeForest.H
├── T8codeForest.C
├── ForestToDirectTopoChange.H
├── ForestToDirectTopoChange.C
├── Make/
│   ├── files
│   └── options
└── dynamicMeshDict
```

### Modificados

| Arquivo | Modificação |
| --- | --- |
| `src/dynamicMesh/Allwmake` | Adicionar subdiretório |
| `wmake/rules/General/T8code` | `T8CODE_INCLUDE`, `T8CODE_LIB` |

## 11. Roadmap

| Fase | Conteúdo | Validação |
| --- | --- | --- |
| 1 | 2D serial, quad, callbacks adapt+interpolate | `icoFoam` vs `dynamicRefineFvMesh` |
| 2 | Mapeamento conservativo (linear + volumétrico) | Massa < 0.1% |
| 3 | 3D hex + `interDyMFoam` (VOF) | damBreak vs tutorial |
| 4 | Paralelo MPI (Partition + Ghost) | 4 processos |
| 5 | Elementos arbitrários + Search | gmsh + solidMechanics |

## 12. Riscos

| Risco | Mitigação |
| --- | --- |
| Conversão forest → directTopoChange | Hex-only na Fase 1 |
| Overhead da conversão | Paper: ~15% do runtime |
| SFC ≠ cell ordering OpenFOAM | Manter cellMap consistente |
| polyPatch types desconhecidos | Mapear do blockMesh inicial |
| MPI init conflito | Passar `MPI_COMM_WORLD` do OF |
| Build cmake ↔ wmake | libt8 separado; wmake linka |
