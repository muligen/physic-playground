using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class drawFluid2d : Image
{
    public int numX = 100;
    public int numY = 60;
    public float[] smoke_c;
    protected override void OnPopulateMesh(VertexHelper vh)
    {
        Vector2 corner1 = Vector2.zero;
        Vector2 corner2 = Vector2.zero;

        corner1.x = 0f;
        corner1.y = 0f;
        corner2.x = 1f;
        corner2.y = 1f;

        corner1.x -= rectTransform.pivot.x;
        corner1.y -= rectTransform.pivot.y;
        corner2.x -= rectTransform.pivot.x;
        corner2.y -= rectTransform.pivot.y;

        corner1.x *= rectTransform.rect.width;
        corner1.y *= rectTransform.rect.height;
        corner2.x *= rectTransform.rect.width;
        corner2.y *= rectTransform.rect.height;

        var uv = Vector4.zero;

        vh.Clear();

        if(smoke_c == null)
            smoke_c = new float[numX * numY];

        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                var c = new Color(smoke_c[i * numY + j], smoke_c[i * numY + j], smoke_c[i * numY + j]);
                vh.AddVert(new Vector3(i * 10, j * 10), c, new Vector2(uv.x, uv.y));
            }
        }
        for (int i = 0; i < numX - 1; i++)
        {
            for (int j = 0; j < numY - 1; j++)
            {
                int idx = i * numY + j;
                vh.AddTriangle(idx, idx + numY, idx + numY+1);
                vh.AddTriangle(idx, idx + numY+1, idx + 1);
            }
        }

    }
}
