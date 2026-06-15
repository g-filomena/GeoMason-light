# GeoMason-light Javadoc style guide

Use this guide for public and protected classes, constructors, fields, and methods.

## Required structure

Every public class should have a short first sentence followed by one or more HTML paragraphs:

```java
/**
 * Represents a directed or undirected spatial edge in a graph.
 *
 * <p>The edge stores the source geometry, endpoint nodes, optional attributes, and cached metrics
 * used by routing and graph-analysis workflows.</p>
 *
 * <h2>Thread-safety</h2>
 * <p>This class is mutable and is not thread-safe.</p>
 *
 * @see Graph
 * @see NodeGraph
 */
```

## HTML conventions

Prefer real Javadoc HTML instead of markdown inside comments:

- Use `<p>...</p>` for paragraphs.
- Use `<ul>` and `<li>` for lists.
- Use `<h2>...</h2>` for longer class-level sections.
- Use `{@code value}` for literals and code fragments.
- Use `{@link ClassName}` or `{@link ClassName#methodName(Type)}` for API links.
- Escape raw generics in prose as `{@code Map<String, AttributeValue>}` rather than writing angle brackets directly.

## Method-level rules

Every non-trivial public method should document:

- what the method returns or mutates;
- whether inputs may be `null`;
- whether coordinates must already be in the same CRS;
- whether the method mutates internal graph/layer state;
- any exception that callers can reasonably recover from.

Example:

```java
/**
 * Returns the first attribute value stored under {@code attributeName}.
 *
 * <p>The lookup is case-sensitive. A missing attribute returns {@code null}; it does not create a new
 * entry in the underlying attribute map.</p>
 *
 * @param attributeName name of the attribute to read; must not be {@code null}
 * @return the stored attribute value, or {@code null} when the attribute is absent
 */
```

## Package-level documentation

Each package should have a `package-info.java` file with:

- a package-level description;
- a short list of central classes;
- notes about CRS, mutability, or graph/routing assumptions;
- `@see` links to adjacent packages.

The package docs in this bundle follow that pattern and can be used as the baseline.
