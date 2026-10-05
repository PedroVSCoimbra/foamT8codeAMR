# Fase 1 — t8code: hexa serial

## Objetivo

Substituir o `hexRef8` pelo t8code em **serial** com
**hexaedros**, mantendo o resto da infraestrutura do
foam-extend intacta. Sem MPI, sem load balancing, sem ghost
layer.

## Arquitetura (minimal)

```text
dynamicFvMesh
  └── dynamicT8codeFvMesh (serial, hex-only)
```

Diferente do documento completo, aqui herdamos diretamente de
`dynamicFvMesh` (não de `topoChangerFvMesh`) — não precisamos
de load balancing na Fase 1.

### Estratégia

O t8code decide **quais** células refinar/desrefinar. Os
comandos são passados para o `directTopoChange`, que produz o
`mapPolyMesh` normalmente. O `fvMesh::updateMesh()` mapeia os
campos automaticamente.

```text
t8code::forest (decide refino)
  → directTopoChange::changeMesh() (executa split/merge)
    → mapPolyMesh
      → fvMesh::updateMesh() (mapeia campos + topologia)
        → correctFluxes (phi = U·Sf nas faces novas)
```

## Fluxo detalhado

```text
mesh.update() ← dynamicT8codeFvMesh
  │
  ├── 1. Converter malha OpenFOAM → forest t8code
  │     (apenas no primeiro passo)
  │
  ├── 2. Se é tempo de refinar (refineInterval):
  │     ├── Copiar campos OF → arrays t8code
  │     │   (alpha1, U, pd como arrays SFC-ordered)
  │     │   t8_forest_set_user_data(forest, &app_data)
  │     │
  │     ├── Adapt: t8_forest_set_adapt() + commit()
  │     │   Callback: refina onde 0.001 < alpha1 < 0.999
  │     │
  │     ├── Balance 2:1: t8_forest_set_balance() + commit()
  │     │
  │     ├── Interpolate: migrate campos entre as forests
  │     │   (injeção simples na Fase 1)
  │     │
  │     └── Converter forest → directTopoChange
  │         Para cada elemento local:
  │           → polyAddPoint / polyAddFace / polyAddCell
  │           → hexRef8::setRefinement era assim
  │
  ├── 3. directTopoChange::changeMesh() + updateMesh()
  │
  └── 4. correctFluxes (phi = U · Sf nas faces novas)
```

## Callbacks

### Adapt

```c
int adapt_cb(forest, forest_from, which_tree, tree_class,
             lelement_id, scheme, is_family, num_elem, elements[])
{
    AppData *d = (AppData*) t8_forest_get_user_data(forest);
    double a1 = d->alpha1[lelement_id];

    if (is_family && a1 < 10.0 && num_elem == 8)
        return -1;                       // coarsen
    if (a1 > 0.001 && a1 < 0.999)
        return 1;                        // refine
    return 0;                            // keep
}
```

### Interpolate (replace)

```c
void replace_cb(forest_old, forest_new, which_tree,
                num_outgoing, first_outgoing,
                num_incoming, first_incoming, user_data)
{
    AppData *d = (AppData*) user_data;

    if (num_outgoing == 1 && num_incoming == 8) {
        // Refino 1→8: injeção simples (Fase 1)
        double val = d->alpha1[first_outgoing[0]];
        for (int i = 0; i < 8; i++)
            d->alpha1[first_incoming[i]] = val;
    }
    else if (num_outgoing == 8 && num_incoming == 1) {
        // Coarsen 8→1: média uniforme (Fase 1)
        double sum = 0;
        for (int i = 0; i < 8; i++)
            sum += d->alpha1[first_outgoing[i]];
        d->alpha1[first_incoming[0]] = sum / 8.0;
    }
}
```

## API t8code usada (serial, sem ghost/partition)

| Função | Para quê |
| --- | --- |
| `t8_forest_init()` | Criar forest |
| `t8_forest_set_cmesh()` | Definir coarse mesh |
| `t8_forest_set_scheme()` | Escolher tipo de elemento |
| `t8_forest_set_level()` | Nível uniforme inicial |
| `t8_forest_set_user_data()` | Passar arrays de campos |
| `t8_forest_set_adapt()` | Callback de refino |
| `t8_forest_set_balance()` | Balanceamento 2:1 |
| `t8_forest_commit()` | Finalizar forest |
| `t8_forest_get_local_num_leaf_elements()` | Nº de elementos |
| `t8_forest_get_leaf_element()` | Obter elemento por ID |
| `t8_forest_element_coordinate()` | Coordenadas do vértice |
| `t8_forest_element_face_neighbor()` | Vizinho de face |
| `t8_forest_element_volume()` | Volume do elemento |
| `t8_element_init()` | Inicializar elemento |
| `t8_element_get_num_faces()` | Nº de faces |
| `t8_element_get_face_corner()` | Vértices de uma face |
| `t8_forest_iterate_replace()` | Migrar dados entre forests |

## Conversão coarse mesh → t8code

```text
polyMesh (blockMesh)
  │
  ├── blockMesh só produz hexaedros
  ├── Para cada célula: verificar se é hex
  │
  ├── Criar t8_cmesh com T8_ECLASS_HEX
  ├── Mapear vértices (coords) e conectividade
  │
  └── t8_forest_set_cmesh()
```

## Conversão forest → directTopoChange

```text
Para cada elemento na forest:
  1. Vértices: t8_forest_element_coordinate()
     → polyAddPoint (com hash merging de pontos partilhados)

  2. Faces: t8_element_get_num_faces()
     + t8_element_get_face_corner()
     → polyAddFace com owner/neighbour

  3. Conexões: t8_forest_element_face_neighbor()
     → owner/neighbour
     → -1 no neighbor = patch de contorno
```

A saída alimenta `directTopoChange::changeMesh()`, que produz
o `mapPolyMesh` como o `hexRef8` fazia.

## Dicionário

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
}
```

## Arquivos

```text
src/dynamicMesh/dynamicFvMesh/dynamicT8codeFvMesh/
├── dynamicT8codeFvMesh.H
├── dynamicT8codeFvMesh.C (update, refine, unrefine)
├── T8codeInitialiser.H  (polyMesh → cmesh)
├── T8codeInitialiser.C
├── T8codeAdaptCallback.H
├── T8codeAdaptCallback.C
├── T8codeReplaceCallback.H
├── T8codeReplaceCallback.C
├── ForestToPolyMesh.H    (forest → directTopoChange)
├── ForestToPolyMesh.C
├── Make/
│   ├── files
│   └── options
└── dynamicMeshDict
```

## Validação

Depois de implementar, rodar o caso `damBreak` do foam-extend
com `dynamicT8codeFvMesh` e comparar com `dynamicRefineFvMesh`:

- Número de células em cada passo
- Posição da frente de onda
- Conservação de massa (alpha1)
- Tempo de CPU por passo

## Dependência do directTopoChange

O `hexRef8::setRefinement()` preenche o `directTopoChange` com
comandos:

```cpp
meshMod.addPoint(/*coords*/);     // ponto novo (meio da aresta)
meshMod.addFace(/*vertices*/);    // face nova (split)
meshMod.addCell(/*master*/);      // célula filha
```

A nossa conversão forest → `directTopoChange` faz o mesmo —
apenas a fonte dos comandos muda (t8code em vez de hexRef8).

```cpp
// Antes (hexRef8):
meshCutter_.setRefinement(cellsToRefine, meshMod);

// Depois (t8code):
T8code::setDirectTopoChange(forest, meshMod);
// interna: para cada elemento novo, polyAddPoint/AddFace/AddCell
directTopoChange::changeMesh(*this, false);
```
