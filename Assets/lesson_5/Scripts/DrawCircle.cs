using JetBrains.Annotations;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class DrawCircle : Image
{
    public List<Vector2> GetCirclePoints(int cnt, float x0, float y0, float x1, float y1)
    {
        var points = new List<Vector2>();
        var mid = new Vector2((x1 + x0) / 2, (y1 + y0) / 2);
        var r = Mathf.Min((x1 - x0) / 2, (y1 - y0)/2);
        var dw = 2 * Mathf.PI / cnt;
        points.Add(mid);
        for(int i = 0; i < cnt; i++)
        {
            points.Add(new Vector2(mid.x + r * Mathf.Cos(i * dw), mid.y + r * Mathf.Sin(i * dw)));
        }
        return points;
    }
    public List<Vector2> GetCirclePointsWithoutCenter(int cnt, float x0, float y0, float x1, float y1)
    {
        var points = new List<Vector2>();
        var mid = new Vector2((x1 + x0) / 2, (y1 + y0) / 2);
        var r = Mathf.Min((x1 - x0) / 2, (y1 - y0)/2);
        var dw = 2 * Mathf.PI / cnt;
        for(int i = 0; i < cnt; i++)
        {
            points.Add(new Vector2(mid.x + r * Mathf.Cos(i * dw), mid.y + r * Mathf.Sin(i * dw)));
        }
        return points;
    }
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

        var points1 = GetCirclePointsWithoutCenter(50, corner1.x, corner1.y, corner2.x, corner2.y);
        float thickness = 10f;
        var points2 = GetCirclePointsWithoutCenter(50, corner1.x + thickness, corner1.y + thickness, corner2.x - thickness, corner2.y - thickness);

        var uv = Vector4.zero;

        var color32 = color;
        vh.Clear();
        //顺时针绘制圆形
        //vh.AddVert(new Vector3(points[0].x, points[0].y), color32, new Vector2(uv.x, uv.y));
        //vh.AddVert(new Vector3(points[1].x, points[1].y), color32, new Vector2(uv.x, uv.y));
        //Debug.Log(points.Count);
        //for (int i = 2; i < points.Count; i++)
        //{
        //    vh.AddVert(new Vector3(points[i].x, points[i].y), color32, new Vector2(uv.x, uv.y));
        //    vh.AddTriangle(i, i-1, 0);
        //}
        //vh.AddTriangle(1, points.Count - 1, 0);

        //顺时针绘制圆环
        vh.AddVert(new Vector3(points1[0].x, points1[0].y), color32, new Vector2(uv.x, uv.y));
        vh.AddVert(new Vector3(points2[0].x, points2[0].y), color32, new Vector2(uv.x, uv.y));
        for (int i = 1; i < points1.Count; i++)
        {
            vh.AddVert(new Vector3(points1[i].x, points1[i].y), color32, new Vector2(uv.x, uv.y));
            vh.AddVert(new Vector3(points2[i].x, points2[i].y), color32, new Vector2(uv.x, uv.y));
            int max_idx = 2 * i + 1;
            vh.AddTriangle(max_idx, max_idx - 1, max_idx - 2);
            vh.AddTriangle(max_idx - 3, max_idx - 2, max_idx - 1);
        }
        int idx = 2 * points1.Count - 1;
        vh.AddTriangle(idx, 1, 0);
        vh.AddTriangle(idx - 1, idx, 0);
    }



}
